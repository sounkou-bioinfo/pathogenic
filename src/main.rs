use clap::Parser;
use csv;
use flate2::read::MultiGzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use reqwest::blocking;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fmt;
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use num_cpus;
use chrono::prelude::*;

/// Command-line arguments
#[derive(Parser)]
#[command(
    name = "pathogenic",
    about = "Identify pathogenic variants from a VCF using official ClinVar (GRCh37 or GRCh38)"
)]
struct Args {
    /// Genome build: must be "GRCh37" (hg19) or "GRCh38" (hg38)
    #[arg(short, long)]
    build: String,

    /// Path to the input VCF (uncompressed)
    #[arg(short, long, value_name = "FILE")]
    input: PathBuf,
}

/// Custom error type for downloads
#[derive(Debug)]
enum DownloadError {
    Io(std::io::Error),
    Reqwest(reqwest::Error),
}
impl fmt::Display for DownloadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DownloadError::Io(e) => write!(f, "IO error: {e}"),
            DownloadError::Reqwest(e) => write!(f, "Reqwest error: {e}"),
        }
    }
}
impl Error for DownloadError {}
impl From<std::io::Error> for DownloadError {
    fn from(e: std::io::Error) -> Self {
        DownloadError::Io(e)
    }
}
impl From<reqwest::Error> for DownloadError {
    fn from(e: reqwest::Error) -> Self {
        DownloadError::Reqwest(e)
    }
}

/// Download a remote file with a progress bar, if not already present locally.
fn download_file(url: &str, out_path: &Path, log_file: &mut File) -> Result<(), DownloadError> {
    println!("  -> Downloading from {url} ...");
    writeln!(log_file, "  -> Downloading from {url} ...")?;
    let mut response = blocking::get(url)?;

    let total_size = response
        .headers()
        .get(reqwest::header::CONTENT_LENGTH)
        .and_then(|ct_len| ct_len.to_str().ok())
        .and_then(|ct_len_str| ct_len_str.parse().ok())
        .unwrap_or(0);

    let pb = ProgressBar::new(total_size);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})",
            )
            .unwrap()
            .progress_chars("=>-"),
    );

    let mut file = File::create(out_path)?;
    let mut downloaded: u64 = 0;
    let mut buffer = [0u8; 8192];
    loop {
        let bytes_read = response.read(&mut buffer)?;
        if bytes_read == 0 {
            break;
        }
        file.write_all(&buffer[..bytes_read])?;
        downloaded += bytes_read as u64;
        pb.set_position(downloaded);
    }
    pb.finish_with_message("Download complete");

    Ok(())
}

/// A specialized type for validated ClinVar lines, storing relevant data
#[derive(Debug, Clone)]
struct ClinVarRecord {
    chr: String,
    pos: u32,
    ref_allele: String,
    alt_allele: String,
    clnsig: String,
    // True if ALT is pathogenic, false if REF is pathogenic and ALT is benign/protective
    is_alt_pathogenic: bool,
    gene: Option<String>,
    allele_id: Option<i32>,
    clnrevstat: Option<String>,
    af_esp: Option<f64>,
    af_exac: Option<f64>,
    af_tgp: Option<f64>,
}

/// Container for ClinVar variants keyed by (chr, pos, ref, alt)
type ClinVarMap = HashMap<(String, u32, String, String), ClinVarRecord>;

/// Parse the semicolon-delimited INFO field into a HashMap
fn parse_info_field(info_str: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for item in info_str.split(';') {
        if item.is_empty() {
            continue;
        }
        let mut kv = item.splitn(2, '=');
        let key = kv.next().unwrap_or("").to_string();
        let val = kv.next().unwrap_or("").to_string();
        map.insert(key, val);
    }
    map
}

/// Convert ClinVar review status to a star rating
fn review_status_to_stars(status: Option<&str>) -> u8 {
    match status {
        Some(s) if s.contains("practice guideline") => 4,
        Some(s) if s.contains("reviewed by expert panel") => 3,
        Some(s) if s.contains("criteria provided, multiple submitters, no conflicts") => 2,
        Some(s) if s.contains("criteria provided, single submitter") => 1,
        _ => 0,
    }
}

/// Parse a ClinVar VCF line into zero or more ClinVarRecords
fn parse_clinvar_line(
    line: &str,
    clinvar_has_chr: bool,
    input_uses_chr: bool,
) -> Option<Vec<ClinVarRecord>> {
    if line.starts_with('#') || line.trim().is_empty() {
        return None;
    }

    let mut fields = line.split('\t');
    let chrom = fields.next()?.to_string();
    let pos_str = fields.next()?;
    let _ = fields.next(); // ID
    let ref_allele = fields.next()?;
    let alt_allele = fields.next()?;
    let _ = fields.next(); // QUAL
    let _ = fields.next(); // FILTER
    let info_str = fields.next()?;
    let pos_num: u32 = pos_str.parse().ok()?;

    let mut chr_fixed = chrom;
    if input_uses_chr && !clinvar_has_chr {
        chr_fixed = if chr_fixed == "MT" {
            "chrM".to_string()
        } else {
            format!("chr{}", chr_fixed)
        };
    } else if !input_uses_chr && clinvar_has_chr {
        chr_fixed = if chr_fixed.eq_ignore_ascii_case("chrM") || chr_fixed.eq_ignore_ascii_case("chrMT") {
            "MT".to_string()
        } else if let Some(stripped) = chr_fixed.strip_prefix("chr") {
            stripped.to_string()
        } else {
            chr_fixed
        };
    }

    let alt_list: Vec<&str> = alt_allele.split(',').collect();
    let info_map = parse_info_field(info_str);

    let clnsig_opt = info_map.get("CLNSIG");
    if clnsig_opt.is_none() {
        return None;
    }
    let clnsig_str = clnsig_opt.unwrap();
    // Only include variants where REF-to-ALT change is pathogenic or likely pathogenic
    if !clnsig_str.contains("Pathogenic") && !clnsig_str.contains("Likely_pathogenic") {
        return None;
    }
    // Exclude ambiguous or conflicting classifications
    if clnsig_str.contains("Conflicting_interpretations_of_pathogenicity") || 
       clnsig_str.contains("Uncertain_significance") {
        return None;
    }
    // Store whether ALT is the pathogenic state (true) or REF is (false)
    let is_alt_pathogenic = !(clnsig_str.contains("Benign") || clnsig_str.contains("Protective"));

    let gene_opt = info_map.get("GENEINFO").map(|g| g.split(':').next().unwrap_or(g).to_string());
    let allele_id_opt = info_map.get("ALLELEID").and_then(|a| a.parse::<i32>().ok());
    let clnrevstat = info_map.get("CLNREVSTAT").map(|s| s.to_string());
    let af_esp = info_map.get("AF_ESP").and_then(|s| s.parse::<f64>().ok());
    let af_exac = info_map.get("AF_EXAC").and_then(|s| s.parse::<f64>().ok());
    let af_tgp = info_map.get("AF_TGP").and_then(|s| s.parse::<f64>().ok());

    let mut recs = Vec::with_capacity(alt_list.len());
    for alt_a in alt_list {
        recs.push(ClinVarRecord {
            chr: chr_fixed.clone(),
            pos: pos_num,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_a.to_string(),
            clnsig: clnsig_str.to_string(),
            is_alt_pathogenic,
            gene: gene_opt.clone(),
            allele_id: allele_id_opt,
            clnrevstat: clnrevstat.clone(),
            af_esp,
            af_exac,
            af_tgp,
        });
    }
    Some(recs)
}

/// Parse the ClinVar .vcf.gz file in parallel
fn parse_clinvar_vcf_gz(
    path_gz: &Path,
    input_uses_chr: bool,
    log_file: &mut File,
) -> Result<(ClinVarMap, String), Box<dyn Error>> {
    println!("\n[STEP] Parsing ClinVar .vcf.gz: {}", path_gz.display());
    writeln!(log_file, "\n[STEP] Parsing ClinVar .vcf.gz: {}", path_gz.display())?;

    let f = File::open(path_gz)?;
    let decoder = MultiGzDecoder::new(f);
    let reader = BufReader::new(decoder);
    // Process lines in parallel without full decompression into memory
    let lines: Vec<String> = reader.lines().collect::<Result<_, _>>()?;
    let total = lines.len() as u64;

    println!(
        "  -> Streaming ClinVar .vcf.gz, processing {} lines with {} cores",
        total,
        num_cpus::get()
    );
    writeln!(
        log_file,
        "  -> Streaming ClinVar .vcf.gz, processing {} lines with {} cores",
        total,
        num_cpus::get()
    )?;

    let mut clinvar_file_date = String::new();
    let mut clinvar_has_chr = false;
    for line in &lines {
        if line.starts_with('#') {
            if line.starts_with("##fileDate") {
                clinvar_file_date = line.trim_start_matches("##fileDate=").to_string();
                break;
            }
            continue;
        }
        let chrom_part = line.split('\t').next().unwrap_or("");
        if chrom_part.starts_with("chr") {
            clinvar_has_chr = true;
        }
        break;
    }
    println!("  -> ClinVar file date: {}", clinvar_file_date);
    writeln!(log_file, "  -> ClinVar file date: {}", clinvar_file_date)?;

    let pb = ProgressBar::new(total);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=>-"),
    );

    let chunk_maps: Vec<ClinVarMap> = lines
        .par_iter()
        .map(|line| {
            pb.inc(1);
            match parse_clinvar_line(line, clinvar_has_chr, input_uses_chr) {
                None => HashMap::new(),
                Some(records) => {
                    let mut local_map = HashMap::with_capacity(records.len());
                    for r in records {
                        let key = (r.chr.clone(), r.pos, r.ref_allele.clone(), r.alt_allele.clone());
                        local_map.insert(key, r);
                    }
                    local_map
                }
            }
        })
        .collect();

    pb.finish_with_message("ClinVar parse complete.");

    let mut final_map = HashMap::new();
    for cm in chunk_maps {
        final_map.extend(cm);
    }

    println!("  -> Final ClinVar map size: {}", final_map.len());
    writeln!(log_file, "  -> Final ClinVar map size: {}", final_map.len())?;
    Ok((final_map, clinvar_file_date))
}

/// Struct for user input variants
#[derive(Debug)]
struct InputVariant {
    chr: String,
    pos: u32,
    ref_allele: String,
    alts: Vec<(String, bool)>, // (alt_allele, is_present_in_genotype)
    genotype: String,
}

/// Parse a user input VCF line
fn parse_input_line(line: &str, user_has_chr: bool, need_chr: bool) -> Option<(String, InputVariant)> {
    if line.starts_with('#') || line.trim().is_empty() {
        return None;
    }

    let mut fields = line.split('\t');
    let chrom = fields.next()?.to_string();
    let pos_str = fields.next()?;
    let _ = fields.next(); // ID
    let ref_allele = fields.next()?;
    let alt_allele = fields.next()?;
    let _ = fields.next(); // QUAL
    let _ = fields.next(); // FILTER
    let _ = fields.next(); // INFO

    let mut rest_cols = Vec::new();
    for c in fields {
        rest_cols.push(c);
    }

    let pos_num: u32 = pos_str.parse().ok()?;
    let mut chr_fixed = chrom;
    if need_chr && !user_has_chr {
        chr_fixed = if chr_fixed == "MT" {
            "chrM".to_string()
        } else {
            format!("chr{}", chr_fixed)
        };
    } else if !need_chr && user_has_chr {
        chr_fixed = if chr_fixed.eq_ignore_ascii_case("chrM") || chr_fixed.eq_ignore_ascii_case("chrMT") {
            "MT".to_string()
        } else if let Some(stripped) = chr_fixed.strip_prefix("chr") {
            stripped.to_string()
        } else {
            chr_fixed
        };
    }

    let alt_list: Vec<String> = alt_allele.split(',').map(|s| s.to_string()).collect();
    let mut present_flags = HashSet::new();

    let genotype = if !rest_cols.is_empty() {
        let format_str = rest_cols[0];
        let format_items: Vec<&str> = format_str.split(':').collect();
        if let Some(gt_index) = format_items.iter().position(|&f| f == "GT") {
            if rest_cols.len() > 1 {
                let sample_str = rest_cols[1];
                let sample_items: Vec<&str> = sample_str.split(':').collect();
                if gt_index < sample_items.len() {
                    sample_items[gt_index].to_string()
                } else {
                    "1/1".to_string() // Default homozygous alternate
                }
            } else {
                "1/1".to_string()
            }
        } else {
            "1/1".to_string()
        }
    } else {
        "1/1".to_string()
    };

    if !genotype.is_empty() {
        let split_gt: Vec<&str> = genotype.split(&['/', '|'][..]).collect();
        for g in split_gt {
            if let Ok(idx) = g.parse::<usize>() {
                if idx >= 1 {
                    present_flags.insert(idx);
                }
            }
        }
    }

    let mut alts = Vec::with_capacity(alt_list.len());
    for (i, alt_a) in alt_list.into_iter().enumerate() {
        let alt_idx = i + 1;
        let is_present = present_flags.contains(&alt_idx);
        alts.push((alt_a, is_present));
    }

    Some((
        line.to_string(),
        InputVariant {
            chr: chr_fixed,
            pos: pos_num,
            ref_allele: ref_allele.to_string(),
            alts,
            genotype,
        },
    ))
}

/// Parse the user input VCF (uncompressed)
fn parse_input_vcf(path: &Path, need_chr: bool, log_file: &mut File) -> Result<Vec<(String, InputVariant)>, Box<dyn Error>> {
    println!("\n[STEP] Parsing user input VCF (uncompressed): {}", path.display());
    writeln!(log_file, "\n[STEP] Parsing user input VCF (uncompressed): {}", path.display())?;

    if let Some(ext) = path.extension() {
        if ext == "gz" {
            return Err(format!("User input VCF cannot be compressed: {}", path.display()).into());
        }
    }

    let mut file = File::open(path)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;

    println!(
        "  -> Loaded user input VCF ({} bytes). Using {} CPU cores.",
        contents.len(),
        num_cpus::get()
    );
    writeln!(
        log_file,
        "  -> Loaded user input VCF ({} bytes). Using {} CPU cores.",
        contents.len(),
        num_cpus::get()
    )?;

    let lines: Vec<&str> = contents.lines().collect();
    let total = lines.len() as u64;

    let mut user_has_chr = false;
    for &l in &lines {
        if l.starts_with('#') || l.trim().is_empty() {
            continue;
        }
        let chrom_part = l.split('\t').next().unwrap_or("");
        if chrom_part.starts_with("chr") {
            user_has_chr = true;
        }
        break;
    }

    let pb = ProgressBar::new(total);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=>-"),
    );

    let chunk_variants: Vec<Vec<(String, InputVariant)>> = lines
        .par_iter()
        .map(|&line| {
            pb.inc(1);
            match parse_input_line(line, user_has_chr, need_chr) {
                None => vec![],
                Some((line, iv)) => vec![(line, iv)],
            }
        })
        .collect();

    pb.finish_with_message("User input parse complete.");

    let final_list: Vec<(String, InputVariant)> = chunk_variants.into_iter().flatten().collect();
    println!("  -> Parsed {} variants from user input.", final_list.len());
    writeln!(log_file, "  -> Parsed {} variants from user input.", final_list.len())?;
    Ok(final_list)
}

/// Main function
fn main() -> Result<(), Box<dyn Error>> {
    let mut log_file = OpenOptions::new()
        .create(true)
        .append(true)
        .open("pathogenic.log")?;

    println!("=== Pathogenic Variant Finder (ClinVar .vcf.gz) ===");
    writeln!(log_file, "=== Pathogenic Variant Finder (ClinVar .vcf.gz) ===")?;
    let now: DateTime<Utc> = Utc::now();
    println!("[LOG] Timestamp: {}", now.to_rfc3339());
    writeln!(log_file, "[LOG] Timestamp: {}", now.to_rfc3339())?;

    let args = Args::parse();
    let build = args.build.to_uppercase();
    let input_path = args.input.clone();

    println!("[LOG] Genome Build: {}", build);
    writeln!(log_file, "[LOG] Genome Build: {}", build)?;
    println!("[LOG] Input File: {}", input_path.display());
    writeln!(log_file, "[LOG] Input File: {}", input_path.display())?;

    println!("[STEP] Checking arguments...");
    writeln!(log_file, "[STEP] Checking arguments...")?;
    if build != "GRCH37" && build != "GRCH38" {
        return Err(format!("Genome build must be 'GRCh37' or 'GRCh38'").into());
    }
    if !input_path.exists() {
        return Err(format!("Input VCF not found: {}", input_path.display()).into());
    }

    let (clinvar_url, tbi_url) = match build.as_str() {
        "GRCH37" => (
            "https://ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz",
            "https://ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi",
        ),
        "GRCH38" => (
            "https://ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
            "https://ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi",
        ),
        _ => unreachable!(),
    };

    let clinvar_dir = Path::new("clinvar_data");
    let gz_path = clinvar_dir.join(format!("clinvar_{}.vcf.gz", build));
    let tbi_path = clinvar_dir.join(format!("clinvar_{}.vcf.gz.tbi", build));

    println!("[STEP] Checking local ClinVar data in {}...", clinvar_dir.display());
    writeln!(log_file, "[STEP] Checking local ClinVar data in {}...", clinvar_dir.display())?;
    fs::create_dir_all(&clinvar_dir)?;

    if !gz_path.exists() {
        println!("  -> No local ClinVar gz found. Downloading...");
        writeln!(log_file, "  -> No local ClinVar gz found. Downloading...")?;
        download_file(clinvar_url, &gz_path, &mut log_file)?;
    } else {
        println!("  -> Found local {}", gz_path.display());
        writeln!(log_file, "  -> Found local {}", gz_path.display())?;
        let check_file = File::open(&gz_path)?;
        let mut check_decoder = MultiGzDecoder::new(check_file);
        let mut buffer = [0u8; 1024];
        if let Err(err) = check_decoder.read(&mut buffer) {
            println!("  -> Local ClinVar gz appears corrupt ({err}). Removing and re-downloading...");
            writeln!(log_file, "  -> Local ClinVar gz appears corrupt ({err}). Removing and re-downloading...")?;
            fs::remove_file(&gz_path)?;
            download_file(clinvar_url, &gz_path, &mut log_file)?;
        }
    }
    if !tbi_path.exists() {
        println!("  -> Missing ClinVar index (.tbi). Downloading...");
        writeln!(log_file, "  -> Missing ClinVar index (.tbi). Downloading...")?;
        download_file(tbi_url, &tbi_path, &mut log_file)?;
    } else {
        println!("  -> Found local {}", tbi_path.display());
        writeln!(log_file, "  -> Found local {}", tbi_path.display())?;
    }

    println!("[STEP] Detecting if user input has 'chr' or not...");
    writeln!(log_file, "[STEP] Detecting if user input has 'chr' or not...")?;
    let file_check = File::open(&input_path)?;
    let mut reader = BufReader::new(file_check);
    let mut line_buf = String::new();
    let mut user_has_chr = false;
    while reader.read_line(&mut line_buf)? > 0 {
        if line_buf.starts_with('#') || line_buf.trim().is_empty() {
            line_buf.clear();
            continue;
        }
        let first_col = line_buf.split('\t').next().unwrap_or("");
        if first_col.starts_with("chr") {
            user_has_chr = true;
        }
        break;
    }
    println!("  -> user input uses chr? {}", user_has_chr);
    writeln!(log_file, "  -> user input uses chr? {}", user_has_chr)?;

    let (clinvar_map, clinvar_file_date) = parse_clinvar_vcf_gz(&gz_path, user_has_chr, &mut log_file)?;
    println!("[LOG] ClinVar File Date: {}", clinvar_file_date);
    writeln!(log_file, "[LOG] ClinVar File Date: {}", clinvar_file_date)?;

    let input_variants = parse_input_vcf(&input_path, user_has_chr, &mut log_file)?;

    println!("\n[STEP] Matching input variants vs. ClinVar...");
    writeln!(log_file, "\n[STEP] Matching input variants vs. ClinVar...")?;
    #[derive(Debug)]
    struct OutputRecord {
        chr: String,
        pos: u32,
        ref_allele: String,
        alt_allele: String,
        clnsig: String,
        // True if ALT is pathogenic, false if REF is pathogenic and ALT is benign/protective
        is_alt_pathogenic: bool,
        gene: Option<String>,
        allele_id: Option<i32>,
        genotype: String,
        review_stars: u8,
        af_esp: Option<f64>,
        af_exac: Option<f64>,
        af_tgp: Option<f64>,
        inheritance: String,
    }

    let pb = ProgressBar::new(input_variants.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=>-"),
    );

    let mut results: Vec<OutputRecord> = input_variants
        .par_iter()
        .flat_map_iter(|(_, iv)| {
            let mut found = Vec::new();
            for (alt_a, is_present) in &iv.alts {
                if !is_present {
                    continue;
                }
                let key = (iv.chr.clone(), iv.pos, iv.ref_allele.clone(), alt_a.clone());
                if let Some(cv) = clinvar_map.get(&key) {
                    // Only include if REF-to-ALT change is pathogenic
                    if !cv.is_alt_pathogenic {
                        continue;
                    }
                    let genotype = iv.genotype.clone();
                    let review_stars = review_status_to_stars(cv.clnrevstat.as_deref());
                    let split_gt: Vec<&str> = genotype.split(&['/', '|'][..]).collect();
                    // Infer inheritance based on genotype
                    let inheritance = if split_gt.len() > 0 && split_gt[0] != "0" && (split_gt.len() == 1 || split_gt[0] != split_gt[1]) {
                        // Haploid or heterozygous cases
                        "Likely Dominant".to_string()
                    } else if split_gt.len() == 2 && split_gt[0] == split_gt[1] && split_gt[0] != "0" {
                        "Dominant or Recessive".to_string()
                    } else {
                        "Unknown".to_string()
                    };
                    found.push(OutputRecord {
                        chr: cv.chr.clone(),
                        pos: cv.pos,
                        ref_allele: cv.ref_allele.clone(),
                        alt_allele: cv.alt_allele.clone(),
                        clnsig: cv.clnsig.clone(),
                        is_alt_pathogenic: cv.is_alt_pathogenic,
                        gene: cv.gene.clone(),
                        allele_id: cv.allele_id,
                        genotype,
                        review_stars,
                        af_esp: cv.af_esp,
                        af_exac: cv.af_exac,
                        af_tgp: cv.af_tgp,
                        inheritance,
                    });
                }
            }
            pb.inc(1);
            found
        })
        .collect();

    pb.finish_with_message("Match complete.");

    // Log matched variants after collection to avoid ? in closure
    for rec in &results {
        println!(
            "Matched variant: chr={}, pos={}, REF={}->ALT={}, CLNSIG={}, Genotype={}",
            rec.chr, rec.pos, rec.ref_allele, rec.alt_allele, rec.clnsig, rec.genotype
        );
        writeln!(
            log_file,
            "Matched variant: chr={}, pos={}, REF={}->ALT={}, CLNSIG={}, Genotype={}",
            rec.chr, rec.pos, rec.ref_allele, rec.alt_allele, rec.clnsig, rec.genotype
        ).unwrap_or(());
    }

    println!("[STEP] Sorting {} matched records...", results.len());
    writeln!(log_file, "[STEP] Sorting {} matched records...", results.len())?;
    results.sort_by(|a, b| {
        let c = a.chr.cmp(&b.chr);
        if c == std::cmp::Ordering::Equal {
            a.pos.cmp(&b.pos)
        } else {
            c
        }
    });

    println!("[STEP] Writing CSV to stdout...\n");
    writeln!(log_file, "[STEP] Writing CSV to stdout...\n")?;
    let mut wtr = csv::Writer::from_writer(std::io::stdout());
    wtr.write_record(&[
        "Chromosome",
        "Position",
        "Ref",
        "Alt",
        "ClinicalSignificance",
        "IsAltPathogenic",
        "Gene",
        "ClinVarAlleleID",
        "Genotype",
        "ReviewStars",
        "AF_ESP",
        "AF_EXAC",
        "AF_TGP",
        "Inheritance",
    ])?;
    for rec in &results {
        wtr.write_record(&[
            &rec.chr,
            &rec.pos.to_string(),
            &rec.ref_allele,
            &rec.alt_allele,
            &rec.clnsig,
            &rec.is_alt_pathogenic.to_string(),
            rec.gene.as_deref().unwrap_or(""),
            &rec.allele_id.map(|id| id.to_string()).unwrap_or_default(),
            &rec.genotype,
            &rec.review_stars.to_string(),
            &rec.af_esp.map(|f| f.to_string()).unwrap_or_default(),
            &rec.af_exac.map(|f| f.to_string()).unwrap_or_default(),
            &rec.af_tgp.map(|f| f.to_string()).unwrap_or_default(),
            &rec.inheritance,
        ])?;
    }
    wtr.flush()?;

    println!("Done. Wrote {} variants to CSV, classified by REF-to-ALT change.\n", results.len());
    writeln!(log_file, "Done. Wrote {} variants to CSV, classified by REF-to-ALT change.\n", results.len())?;
    Ok(())
}
