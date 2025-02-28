use clap::Parser;
use csv;
use flate2::read::MultiGzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use reqwest::blocking;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fmt;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};

use num_cpus;

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
fn download_file(url: &str, out_path: &Path) -> Result<(), DownloadError> {
    println!("  -> Downloading from {url} ...");
    let mut response = blocking::get(url)?;

    // If the server provides a Content-Length, we can show progress
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

/// A specialized type for validated ClinVar lines
/// storing only the data we care about: chr, pos, ref, alt, clnsig, gene, allele_id
#[derive(Debug, Clone)]
struct ClinVarRecord {
    chr: String,
    pos: u32,
    ref_allele: String,
    alt_allele: String,
    clnsig: String,
    gene: Option<String>,
    allele_id: Option<i32>,
}

/// Container for final ClinVar variants keyed by (chr, pos, ref, alt)
type ClinVarMap = HashMap<(String, u32, String, String), ClinVarRecord>;

/// A helper to parse the semicolon-delimited INFO field
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

/// Attempt to parse one line of the ClinVar VCF into zero or more `ClinVarRecord`s.
/// - Return `None` if line is invalid or not Pathogenic.
/// - If multiple ALT alleles, return multiple records.
fn parse_clinvar_line(
    line: &str,
    clinvar_has_chr: bool,
    input_uses_chr: bool,
) -> Option<Vec<ClinVarRecord>> {
    // Skip comments
    if line.starts_with('#') || line.trim().is_empty() {
        return None;
    }

    // Format: CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO ...
    let mut fields = line.split('\t');
    let chrom = fields.next()?.to_string();
    let pos_str = fields.next()?;
    // skip ID
    let _ = fields.next();
    let ref_allele = fields.next()?;
    let alt_allele = fields.next()?;
    // skip QUAL, FILTER
    let _ = fields.next();
    let _ = fields.next();
    let info_str = fields.next()?;
    let pos_num: u32 = pos_str.parse().ok()?;

    // unify naming
    let mut chr_fixed = chrom;
    if input_uses_chr && !clinvar_has_chr {
        // add "chr"
        if chr_fixed == "MT" {
            chr_fixed = "chrM".to_string();
        } else {
            chr_fixed = format!("chr{}", chr_fixed);
        }
    } else if !input_uses_chr && clinvar_has_chr {
        // remove "chr"
        if chr_fixed.eq_ignore_ascii_case("chrM") || chr_fixed.eq_ignore_ascii_case("chrMT") {
            chr_fixed = "MT".to_string();
        } else if let Some(stripped) = chr_fixed.strip_prefix("chr") {
            chr_fixed = stripped.to_string();
        }
    }

    let alt_list: Vec<&str> = alt_allele.split(',').collect();
    let info_map = parse_info_field(info_str);

    let clnsig_opt = info_map.get("CLNSIG");
    if clnsig_opt.is_none() {
        return None;
    }
    let clnsig_str = clnsig_opt.unwrap();
    // Must contain "Pathogenic" or "Likely_pathogenic"
    if !clnsig_str.contains("Pathogenic") {
        return None;
    }
    // Skip if "Conflicting_interpretations_of_pathogenicity"
    if clnsig_str.contains("Conflicting_interpretations_of_pathogenicity") {
        return None;
    }

    let gene_opt = info_map
        .get("GENEINFO")
        .map(|g| g.split(':').next().unwrap_or(g).to_string());
    let allele_id_opt = info_map.get("ALLELEID").and_then(|a| a.parse::<i32>().ok());

    // Return one ClinVarRecord for each ALT
    let mut recs = Vec::with_capacity(alt_list.len());
    for alt_a in alt_list {
        recs.push(ClinVarRecord {
            chr: chr_fixed.clone(),
            pos: pos_num,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_a.to_string(),
            clnsig: clnsig_str.to_string(),
            gene: gene_opt.clone(),
            allele_id: allele_id_opt,
        });
    }
    Some(recs)
}

/// Parse the entire ClinVar .vcf.gz in parallel, building a `ClinVarMap`.
fn parse_clinvar_vcf_gz(
    path_gz: &Path,
    input_uses_chr: bool,
) -> Result<ClinVarMap, Box<dyn Error>> {
    println!("\n[STEP] Parsing ClinVar .vcf.gz: {}", path_gz.display());

    // We can't just read the entire .gz into memory easily. Instead, let's read line by line
    // but in big chunk buffers. Then we can parallelize in a chunked manner if we want
    // a "type-driven" approach. For simplicity, let's read line by line in a single pass,
    // then do parallel "map" on lines if we store them in memory.

    // If the file is extremely large, it might be more memory than you'd like, but ClinVar
    // is ~300-500MB compressed, ~3GB uncompressed. We'll do it carefully.

    // 1) Read + decompress line by line into memory
    // 2) Detect if ClinVar has "chr" from the first data line
    // 3) Parallel parse each line with parse_clinvar_line()
    // 4) Merge into final map

    let f = File::open(path_gz)?;
    let mut decoder = MultiGzDecoder::new(f);

    // We'll store lines in a buffer. Because of memory constraints, if you prefer,
    // you can parse line-by-line in a single thread or in a streaming fashion.
    // For demonstration, we'll do the "load into memory" approach:

    let mut unzipped_contents = String::new();
    decoder.read_to_string(&mut unzipped_contents)?;

    println!(
        "  -> Decompressed size: {} bytes. Using {} CPU cores.",
        unzipped_contents.len(),
        num_cpus::get()
    );

    let lines: Vec<&str> = unzipped_contents.lines().collect();
    let total = lines.len() as u64;

    // Does ClinVar have "chr"? We'll check the first data line
    let mut clinvar_has_chr = false;
    for &line in &lines {
        if line.starts_with('#') {
            continue;
        }
        let chrom_part = line.split('\t').next().unwrap_or("");
        if chrom_part.starts_with("chr") {
            clinvar_has_chr = true;
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

    // Parallel parse
    let chunk_maps: Vec<ClinVarMap> = lines
        .par_iter()
        .map(|&line| {
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

    // Merge chunk maps
    let mut final_map = HashMap::new();
    for cm in chunk_maps {
        final_map.extend(cm);
    }

    println!("  -> Final ClinVar map size: {}", final_map.len());
    Ok(final_map)
}

//-------------------------------------
// Parse the user input (uncompressed)
//-------------------------------------

/// Simple struct for user input variants
#[derive(Debug)]
struct InputVariant {
    chr: String,
    pos: u32,
    ref_allele: String,
    alts: Vec<(String, bool)>, // (alt_allele, is_present_in_genotype)
}

/// Parse user input (uncompressed) VCF line into zero/one `InputVariant`.
fn parse_input_line(line: &str, user_has_chr: bool, need_chr: bool) -> Option<InputVariant> {
    if line.starts_with('#') || line.trim().is_empty() {
        return None;
    }

    // CHROM POS ID REF ALT QUAL FILTER INFO [FORMAT, SAMPLE...]
    let mut fields = line.split('\t');
    let chrom = fields.next()?.to_string();
    let pos_str = fields.next()?;
    // skip ID
    let _ = fields.next();
    let ref_allele = fields.next()?;
    let alt_allele = fields.next()?;
    // skip QUAL, FILTER
    let _ = fields.next();
    let _ = fields.next();
    // skip INFO
    let _ = fields.next();

    let mut rest_cols = Vec::new();
    for c in fields {
        rest_cols.push(c);
    }

    let pos_num: u32 = pos_str.parse().ok()?;
    // unify naming
    let mut chr_fixed = chrom;
    if need_chr && !user_has_chr {
        if chr_fixed == "MT" {
            chr_fixed = "chrM".to_string();
        } else {
            chr_fixed = format!("chr{}", chr_fixed);
        }
    } else if !need_chr && user_has_chr {
        if chr_fixed.eq_ignore_ascii_case("chrM") || chr_fixed.eq_ignore_ascii_case("chrMT") {
            chr_fixed = "MT".to_string();
        } else if let Some(stripped) = chr_fixed.strip_prefix("chr") {
            chr_fixed = stripped.to_string();
        }
    }

    let alt_list: Vec<String> = alt_allele.split(',').map(|s| s.to_string()).collect();
    let mut present_flags = HashSet::new();

    // If we have at least one col, the first is FORMAT
    // if we have second col, that's the first sample
    if !rest_cols.is_empty() {
        let format_str = rest_cols[0];
        let format_items: Vec<&str> = format_str.split(':').collect();
        let gt_index_opt = format_items.iter().position(|&f| f == "GT");

        if gt_index_opt.is_some() && rest_cols.len() > 1 {
            let sample_str = rest_cols[1];
            let sample_items: Vec<&str> = sample_str.split(':').collect();
            let gt_index = gt_index_opt.unwrap();
            if gt_index < sample_items.len() {
                let gt_val = sample_items[gt_index];
                let split_gt: Vec<&str> = gt_val.split(&['/', '|'][..]).collect();
                for g in split_gt {
                    if let Ok(idx) = g.parse::<usize>() {
                        if idx >= 1 {
                            present_flags.insert(idx);
                        }
                    }
                }
            } else {
                // no GT => assume all alt
                for i in 1..=alt_list.len() {
                    present_flags.insert(i);
                }
            }
        } else {
            // no GT => assume all alt
            for i in 1..=alt_list.len() {
                present_flags.insert(i);
            }
        }
    } else {
        // no format => assume all alt
        for i in 1..=alt_list.len() {
            present_flags.insert(i);
        }
    }

    let mut alts = Vec::with_capacity(alt_list.len());
    for (i, alt_a) in alt_list.into_iter().enumerate() {
        let alt_idx = i + 1;
        let is_present = present_flags.contains(&alt_idx);
        alts.push((alt_a, is_present));
    }

    Some(InputVariant {
        chr: chr_fixed,
        pos: pos_num,
        ref_allele: ref_allele.to_string(),
        alts,
    })
}

/// Parse the entire user input VCF (uncompressed).
fn parse_input_vcf(path: &Path, need_chr: bool) -> Result<Vec<InputVariant>, Box<dyn Error>> {
    println!("\n[STEP] Parsing user input VCF (uncompressed): {}", path.display());

    // We'll do a quick check: if path ends in ".gz", we reject
    if let Some(ext) = path.extension() {
        if ext == "gz" {
            return Err(format!("User input VCF cannot be compressed: {}", path.display()).into());
        }
    }

    // Read entire file into memory for parallel processing
    let mut file = File::open(path)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;

    println!(
        "  -> Loaded user input VCF ({} bytes). Using {} CPU cores.",
        contents.len(),
        num_cpus::get()
    );

    let lines: Vec<&str> = contents.lines().collect();
    let total = lines.len() as u64;

    // Detect if user input has "chr" from the first data line
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

    let chunk_variants: Vec<Vec<InputVariant>> = lines
        .par_iter()
        .map(|&line| {
            pb.inc(1);
            match parse_input_line(line, user_has_chr, need_chr) {
                None => vec![],
                Some(iv) => vec![iv],
            }
        })
        .collect();

    pb.finish_with_message("User input parse complete.");

    let final_list: Vec<InputVariant> = chunk_variants.into_iter().flatten().collect();
    println!("  -> Parsed {} variants from user input.", final_list.len());
    Ok(final_list)
}

//-------------------------------------
// Main
//-------------------------------------

fn main() -> Result<(), Box<dyn Error>> {
    println!("=== Pathogenic Variant Finder (ClinVar .vcf.gz) ===\n");

    let args = Args::parse();
    let build = args.build.to_uppercase();
    let input_path = args.input.clone();

    println!("[STEP] Checking arguments...");
    if build != "GRCH37" && build != "GRCH38" {
        return Err(format!("Genome build must be 'GRCh37' or 'GRCh38'").into());
    }
    if !input_path.exists() {
        return Err(format!("Input VCF not found: {}", input_path.display()).into());
    }

    // ClinVar file URLs
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

    // Local paths
    let clinvar_dir = Path::new("clinvar_data");
    let gz_path = clinvar_dir.join(format!("clinvar_{}.vcf.gz", build));
    let tbi_path = clinvar_dir.join(format!("clinvar_{}.vcf.gz.tbi", build));

    println!("[STEP] Checking local ClinVar data in {}...", clinvar_dir.display());
    fs::create_dir_all(&clinvar_dir)?;

    // Download if missing
    if !gz_path.exists() {
        println!("  -> No local ClinVar gz found. Downloading...");
        download_file(clinvar_url, &gz_path)?;
    } else {
        println!("  -> Found local {}", gz_path.display());
        // Verify that the file is valid by decompressing a small portion.
        // If this fails, remove the file and re-download it.
        {
            let check_file = File::open(&gz_path)?;
            let mut check_decoder = MultiGzDecoder::new(check_file);
            let mut buffer = [0u8; 1024];
            if let Err(err) = check_decoder.read(&mut buffer) {
                println!("  -> Local ClinVar gz appears corrupt ({err}). Removing and re-downloading...");
                fs::remove_file(&gz_path)?;
                download_file(clinvar_url, &gz_path)?;
            }
        }
    }
    if !tbi_path.exists() {
        println!("  -> Missing ClinVar index (.tbi). Downloading...");
        download_file(tbi_url, &tbi_path)?;
    } else {
        println!("  -> Found local {}", tbi_path.display());
    }

    // First, detect if user input wants "chr"
    // We'll do a quick check on the first data line
    println!("[STEP] Detecting if user input has 'chr' or not...");
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

    // Parse the ClinVar .vcf.gz in parallel
    let clinvar_map = parse_clinvar_vcf_gz(&gz_path, user_has_chr)?;

    // Parse the user input (uncompressed)
    let input_variants = parse_input_vcf(&input_path, user_has_chr)?;

    // Matching
    println!("\n[STEP] Matching input variants vs. ClinVar...");
    #[derive(Debug)]
    struct OutputRecord {
        chr: String,
        pos: u32,
        ref_allele: String,
        alt_allele: String,
        clnsig: String,
        gene: Option<String>,
        allele_id: Option<i32>,
    }

    let pb = ProgressBar::new(input_variants.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=>-"),
    );

    // Collect in parallel
    let mut results: Vec<OutputRecord> = input_variants
        .par_iter()
        .flat_map_iter(|iv| {
            // For each alt that is present
            let mut found = Vec::new();
            for (alt_a, is_present) in &iv.alts {
                if !is_present {
                    continue;
                }
                let key = (iv.chr.clone(), iv.pos, iv.ref_allele.clone(), alt_a.clone());
                if let Some(cv) = clinvar_map.get(&key) {
                    found.push(OutputRecord {
                        chr: cv.chr.clone(),
                        pos: cv.pos,
                        ref_allele: cv.ref_allele.clone(),
                        alt_allele: cv.alt_allele.clone(),
                        clnsig: cv.clnsig.clone(),
                        gene: cv.gene.clone(),
                        allele_id: cv.allele_id,
                    });
                }
            }
            pb.inc(1);
            found
        })
        .collect();

    pb.finish_with_message("Match complete.");

    println!("[STEP] Sorting {} matched records...", results.len());
    results.sort_by(|a, b| {
        let c = a.chr.cmp(&b.chr);
        if c == std::cmp::Ordering::Equal {
            a.pos.cmp(&b.pos)
        } else {
            c
        }
    });

    println!("[STEP] Writing CSV to stdout...\n");
    let mut wtr = csv::Writer::from_writer(std::io::stdout());
    wtr.write_record(&[
        "Chromosome",
        "Position",
        "Ref",
        "Alt",
        "ClinicalSignificance",
        "Gene",
        "ClinVarAlleleID",
    ])?;
    for rec in &results {
        wtr.write_record(&[
            &rec.chr,
            &rec.pos.to_string(),
            &rec.ref_allele,
            &rec.alt_allele,
            &rec.clnsig,
            rec.gene.as_deref().unwrap_or(""),
            &rec.allele_id.map(|id| id.to_string()).unwrap_or_default(),
        ])?;
    }
    wtr.flush()?;

    println!("Done. {} total pathogenic variants matched.\n", results.len());
    Ok(())
}
