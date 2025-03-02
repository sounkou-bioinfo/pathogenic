use chrono::prelude::*;
use clap::Parser;
use csv;
use flate2::read::MultiGzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use num_cpus;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fmt;
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use noodles::vcf::variant::record::info::field::Value;
use noodles::vcf::io::indexed_reader::IndexedReader;
use noodles::vcf::variant::record::AlternateBases;
use noodles::core::region::Region;
use noodles::csi;

/// Command-line arguments
#[derive(Parser)]
#[command(
    name = "pathogenic",
    about = "Identify pathogenic variants from a VCF using official ClinVar plus 1000 Genomes frequencies (GRCh37 or GRCh38)"
)]
struct Args {
    /// Genome build: must be "GRCh37" (hg19) or "GRCh38" (hg38)
    #[arg(short, long)]
    build: String,

    /// Path to the input VCF (can be uncompressed or gzipped)
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

/// Download a remote file with a progress bar using parallel chunk downloads,
/// if the server supports HTTP range requests. Otherwise, fall back to a single-threaded download.
fn download_file(
    url: &str,
    out_path: &Path,
    log_file: &mut File,
) -> Result<(), DownloadError> {
    use std::sync::{
        atomic::{AtomicU64, Ordering},
        Arc,
    };
    use std::thread;

    println!("  -> Starting download from {url}");
    writeln!(log_file, "  -> Starting download from {url}")?;

    let client = reqwest::blocking::Client::new();
    let head_resp = client.head(url).send()?;
    let total_size = head_resp
        .headers()
        .get(reqwest::header::CONTENT_LENGTH)
        .and_then(|s| s.to_str().ok())
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(0);
    let accept_ranges = head_resp
        .headers()
        .get(reqwest::header::ACCEPT_RANGES)
        .and_then(|s| s.to_str().ok())
        .unwrap_or("");

    if total_size == 0 || accept_ranges != "bytes" {
        println!("  -> Server does not support parallel downloads; falling back to single-threaded download.");
        writeln!(log_file, "  -> Server does not support parallel downloads; falling back to single-threaded download.")?;
        let mut response = client.get(url).send()?;
        let pb = ProgressBar::new(total_size);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})")
                .unwrap()
                .progress_chars("=>-"),
        );
        let mut file = File::create(out_path)?;
        let mut downloaded: u64 = 0;
        let mut buffer = [0u8; 8192];
        loop {
            let bytes_read = response.read(&mut buffer)?;
            if bytes_read == 0 { break; }
            file.write_all(&buffer[..bytes_read])?;
            downloaded += bytes_read as u64;
            pb.set_position(downloaded);
        }
        pb.finish_with_message("Download complete");
        return Ok(());
    }

    let num_chunks = num_cpus::get();
    let chunk_size = total_size / num_chunks as u64;
    let mut ranges = Vec::new();
    for i in 0..num_chunks {
        let start = i as u64 * chunk_size;
        let end = if i == num_chunks - 1 { total_size - 1 } else { (i as u64 + 1) * chunk_size - 1 };
        ranges.push((start, end));
    }

    let progress = Arc::new(AtomicU64::new(0));
    let mut handles = Vec::new();

    for (i, (start, end)) in ranges.into_iter().enumerate() {
        let client = client.clone();
        let url = url.to_string();
        let progress = Arc::clone(&progress);
        let handle = thread::spawn(move || -> Result<(usize, Vec<u8>), DownloadError> {
            let range_header = format!("bytes={}-{}", start, end);
            let mut resp = client.get(&url)
                .header(reqwest::header::RANGE, range_header)
                .send()?;
            let mut buf = Vec::new();
            let mut local_buf = [0u8; 8192];
            loop {
                let n = resp.read(&mut local_buf)?;
                if n == 0 { break; }
                buf.extend_from_slice(&local_buf[..n]);
                progress.fetch_add(n as u64, Ordering::Relaxed);
            }
            Ok((i, buf))
        });
        handles.push(handle);
    }

    let pb = ProgressBar::new(total_size);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{total_bytes} ({eta})")
            .unwrap()
            .progress_chars("=>-"),
    );

    while progress.load(Ordering::Relaxed) < total_size {
        pb.set_position(progress.load(Ordering::Relaxed));
        std::thread::sleep(std::time::Duration::from_millis(100));
    }
    pb.finish_with_message("Download complete");

    let mut chunks: Vec<(usize, Vec<u8>)> = Vec::with_capacity(handles.len());
    for handle in handles {
        let res = handle.join().map_err(|_| {
            DownloadError::Io(std::io::Error::new(std::io::ErrorKind::Other, "Thread join error"))
        })??;
        chunks.push(res);
    }
    chunks.sort_by_key(|(i, _)| *i);
    let mut file = File::create(out_path)?;
    for (_i, data) in chunks {
        file.write_all(&data)?;
    }
    Ok(())
}

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


/// Custom error to unify reading different VCFs (ClinVar or 1000G)
enum VcfReadError {
    Io(std::io::Error),
    Parse(std::num::ParseIntError),
}

impl From<std::io::Error> for VcfReadError {
    fn from(e: std::io::Error) -> Self {
        VcfReadError::Io(e)
    }
}

impl From<std::num::ParseIntError> for VcfReadError {
    fn from(e: std::num::ParseIntError) -> Self {
        VcfReadError::Parse(e)
    }
}

impl fmt::Display for VcfReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VcfReadError::Io(e) => write!(f, "IO error during VCF read: {}", e),
            VcfReadError::Parse(e) => write!(f, "Parse error during VCF read: {}", e),
        }
    }
}

/// A specialized type for records from the ClinVar VCF
#[derive(Debug, Clone)]
struct ClinVarRecord {
    chr: String,
    pos: u32,
    ref_allele: String,
    alt_allele: String,
    clnsig: String,
    is_alt_pathogenic: bool,
    gene: Option<String>,
    allele_id: Option<i32>,
    clnrevstat: Option<String>,
    clndn: Option<String>,
    af_esp: Option<f64>,
    af_exac: Option<f64>,
    af_tgp: Option<f64>,
}

/// Container for ClinVar variants keyed by (chr, pos, ref, alt)
type ClinVarMap = HashMap<(String, u32, String, String), ClinVarRecord>;

/// Parse a single line from the ClinVar VCF
fn parse_clinvar_line(line: &str) -> Option<Vec<ClinVarRecord>> {
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

    let pos_num = pos_str.parse::<u32>().ok()?;

    let chr_norm = if chrom.eq_ignore_ascii_case("chrM") || chrom.eq_ignore_ascii_case("chrMT") {
        "MT".to_string()
    } else if let Some(stripped) = chrom.strip_prefix("chr") {
        stripped.to_string()
    } else {
        chrom.clone()
    };

    let alt_list: Vec<&str> = alt_allele.split(',').collect();
    let info_map = parse_info_field(info_str);

    let clnsig_opt = info_map.get("CLNSIG");
    if clnsig_opt.is_none() {
        return None;
    }
    let clnsig_str = clnsig_opt.unwrap();

    // Only include variants where REF->ALT is pathogenic or likely pathogenic
    if !clnsig_str.contains("Pathogenic") && !clnsig_str.contains("Likely_pathogenic") {
        return None;
    }

    // Exclude ambiguous or conflicting classifications
    if clnsig_str.contains("Conflicting_interpretations_of_pathogenicity")
        || clnsig_str.contains("Uncertain_significance")
    {
        return None;
    }

    // If CLNSIG suggests benign, alt is not pathogenic
    let is_alt_pathogenic = !(clnsig_str.contains("Benign") || clnsig_str.contains("Protective"));

    let gene_opt = info_map
        .get("GENEINFO")
        .map(|g| g.split(':').next().unwrap_or(g).to_string());
    let allele_id_opt = info_map.get("ALLELEID").and_then(|a| a.parse::<i32>().ok());
    let clnrevstat = info_map.get("CLNREVSTAT").map(|s| s.to_string());
    let clndn = info_map.get("CLNDN").cloned();
    let af_esp = info_map.get("AF_ESP").and_then(|s| s.parse::<f64>().ok());
    let af_exac = info_map.get("AF_EXAC").and_then(|s| s.parse::<f64>().ok());
    let af_tgp = info_map.get("AF_TGP").and_then(|s| s.parse::<f64>().ok());

    let mut recs = Vec::with_capacity(alt_list.len());
    for alt_a in alt_list {
        recs.push(ClinVarRecord {
            chr: chr_norm.clone(),
            pos: pos_num,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_a.to_string(),
            clnsig: clnsig_str.to_string(),
            is_alt_pathogenic,
            gene: gene_opt.clone(),
            allele_id: allele_id_opt,
            clnrevstat: clnrevstat.clone(),
            clndn: clndn.clone(),
            af_esp,
            af_exac,
            af_tgp,
        });
    }
    Some(recs)
}

/// Parse the ClinVar VCF in parallel
fn parse_clinvar_vcf_gz(
    path_gz: &Path,
    log_file: &mut File,
) -> Result<(ClinVarMap, String), Box<dyn Error>> {
    println!("[STEP] Parsing ClinVar VCF: {}", path_gz.display());
    writeln!(log_file, "[STEP] Parsing ClinVar VCF: {}", path_gz.display())?;

    let file = File::open(path_gz)?;
    let decoder = MultiGzDecoder::new(file);
    let reader = BufReader::new(decoder);

    let lines: Vec<String> = reader.lines().collect::<Result<_, _>>()?;
    let total = lines.len() as u64;

    println!("  -> Streaming ClinVar VCF, total lines: {}", total);
    writeln!(
        log_file,
        "  -> Streaming ClinVar VCF, total lines: {}",
        total
    )?;

    let mut file_date = String::new();
    for line in &lines {
        if line.starts_with("##fileDate") {
            file_date = line.trim_start_matches("##fileDate=").to_string();
            break;
        }
    }

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
            match parse_clinvar_line(line) {
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
    Ok((final_map, file_date))
}

/// A specialized type for user input variants
#[derive(Debug)]
struct InputVariant {
    chr: String,
    pos: u32,
    ref_allele: String,
    alts: Vec<(String, bool)>, // (alt_allele, is_present_in_genotype)
    genotype: String,
}

/// Attempt to parse a single user line from VCF
fn parse_input_line(line: &str) -> Option<(String, InputVariant)> {
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

    let pos_num = pos_str.parse::<u32>().ok()?;
    let chr_fixed = if chrom.eq_ignore_ascii_case("chrM") || chrom.eq_ignore_ascii_case("chrMT") {
        "MT".to_string()
    } else if let Some(stripped) = chrom.strip_prefix("chr") {
        stripped.to_string()
    } else {
        chrom.to_string()
    };

    let alt_list: Vec<String> = alt_allele.split(',').map(|s| s.to_string()).collect();
    let mut present_flags = HashSet::new();

    // Attempt to parse genotype from the next columns
    // Usually the first of rest_cols is the FORMAT, second is the sample data
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
                    "1/1".to_string() 
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

    Some((line.to_string(), InputVariant {
        chr: chr_fixed,
        pos: pos_num,
        ref_allele: ref_allele.to_string(),
        alts,
        genotype,
    }))
}

/// Parse the user input VCF (can be gzipped or uncompressed)
fn parse_input_vcf(path: &Path, log_file: &mut File) -> Result<Vec<(String, InputVariant)>, Box<dyn Error>> {
    println!(
        "[STEP] Parsing user input VCF: {}",
        path.display()
    );
    writeln!(
        log_file,
        "[STEP] Parsing user input VCF: {}",
        path.display()
    )?;

    if !path.exists() {
        return Err(format!("User input VCF not found: {}", path.display()).into());
    }

    // Check if path is gzipped by extension
    let mut lines = Vec::new();
    if let Some(ext) = path.extension() {
        if ext == "gz" {
            let file = File::open(path)?;
            let decoder = MultiGzDecoder::new(file);
            let reader = BufReader::new(decoder);
            for l in reader.lines() {
                lines.push(l?);
            }
        } else {
            let file = File::open(path)?;
            let reader = BufReader::new(file);
            for l in reader.lines() {
                lines.push(l?);
            }
        }
    } else {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        for l in reader.lines() {
            lines.push(l?);
        }
    }

    println!(
        "  -> Loaded user VCF ({} lines). Using {} CPU cores.",
        lines.len(),
        num_cpus::get()
    );
    writeln!(
        log_file,
        "  -> Loaded user VCF ({} lines). Using {} CPU cores.",
        lines.len(),
        num_cpus::get()
    )?;

    let total = lines.len() as u64;
    let pb = ProgressBar::new(total);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=>-"),
    );

    let chunk_variants: Vec<Vec<(String, InputVariant)>> = lines
        .par_iter()
        .map(|line| {
            pb.inc(1);
            match parse_input_line(line) {
                None => vec![],
                Some((l, iv)) => vec![(l, iv)],
            }
        })
        .collect();

    pb.finish_with_message("User input parse complete.");

    let final_list: Vec<(String, InputVariant)> = chunk_variants.into_iter().flatten().collect();
    println!("  -> Parsed {} variants from user input.", final_list.len());
    writeln!(
        log_file,
        "  -> Parsed {} variants from user input.",
        final_list.len()
    )?;
    Ok(final_list)
}

/// A specialized type for records from the ClinVar pathogenic summary TSV.
/// These fields are used for annotation in the final output.
#[derive(Debug, Clone)]
struct ClinVarTsvRecord {
    molecular_consequence: Option<String>,
    functional_consequence: Option<String>,
    mode_of_inheritance: Option<String>,
    preferred_values: Option<String>,
    citations: Option<String>,
    comments: Option<String>,
    family_data: Option<String>,
    record_status: Option<String>,
    description: Option<String>,
    date_last_evaluated: Option<String>,
}

/// Parse the ClinVar pathogenic summary TSV for deeper annotation
fn parse_clinvar_tsv(
    tsv_path: &Path,
    build: &str,
    log_file: &mut File,
) -> Result<HashMap<(String, u32, String, String), ClinVarTsvRecord>, Box<dyn Error>> {
    println!("[STEP] Parsing ClinVar summary TSV: {}", tsv_path.display());
    writeln!(
        log_file,
        "[STEP] Parsing ClinVar summary TSV: {}",
        tsv_path.display()
    )?;

    let file = File::open(tsv_path)?;
    let reader: Box<dyn BufRead> = if let Some(ext) = tsv_path.extension() {
        if ext == "xz" {
            Box::new(BufReader::new(xz2::read::XzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        }
    } else {
        Box::new(BufReader::new(file))
    };

    let mut line_count = 0;
    let mut map = HashMap::new();
    let mut header_indices: Option<HashMap<String, usize>> = None;

    // We handle GRCh37 or GRCh38 columns
    let want_37 = build == "GRCH37";

    for line_result in reader.lines() {
        let line = line_result?;
        line_count += 1;

        if line_count == 1 {
            // Header row
            let headers: Vec<&str> = line.split('\t').collect();
            let mut idx_map = HashMap::new();
            for (i, h) in headers.iter().enumerate() {
                idx_map.insert(h.to_string(), i);
            }
            header_indices = Some(idx_map);
            continue;
        }
        let headers_map = header_indices.as_ref().ok_or("No header row found")?;
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 28 {
            continue;
        }

        let c_chr = if want_37 {
            headers_map.get("GRCh37_Chromosome")
        } else {
            headers_map.get("GRCh38_Chromosome")
        };
        let c_pos = if want_37 {
            headers_map.get("GRCh37_PositionVCF")
        } else {
            headers_map.get("GRCh38_PositionVCF")
        };
        let c_ref = if want_37 {
            headers_map.get("GRCh37_ReferenceAlleleVCF")
        } else {
            headers_map.get("GRCh38_ReferenceAlleleVCF")
        };
        let c_alt = if want_37 {
            headers_map.get("GRCh37_AlternateAlleleVCF")
        } else {
            headers_map.get("GRCh38_AlternateAlleleVCF")
        };

        if c_chr.is_none() || c_pos.is_none() || c_ref.is_none() || c_alt.is_none() {
            continue;
        }
        let chr_val = cols[*c_chr.unwrap()].trim();
        let pos_val = cols[*c_pos.unwrap()].trim();
        let ref_val = cols[*c_ref.unwrap()].trim();
        let alt_val = cols[*c_alt.unwrap()].trim();

        let pos_parsed = pos_val.parse::<u32>().unwrap_or(0);
        if chr_val.is_empty() || pos_parsed == 0 || ref_val.is_empty() || alt_val.is_empty() {
            continue;
        }

        // Normalize chromosome
        let chr_norm = if chr_val.eq_ignore_ascii_case("chrM")
            || chr_val.eq_ignore_ascii_case("chrMT")
        {
            "MT".to_string()
        } else if let Some(stripped) = chr_val.strip_prefix("chr") {
            stripped.to_string()
        } else {
            chr_val.to_string()
        };

        let rec = ClinVarTsvRecord {
            molecular_consequence: headers_map
                .get("MolecularConsequence")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            functional_consequence: headers_map
                .get("FunctionalConsequence")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            mode_of_inheritance: headers_map
                .get("ModeOfInheritance")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            preferred_values: headers_map
                .get("PreferredValues")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            citations: headers_map
                .get("Citations")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            comments: headers_map
                .get("Comments")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            family_data: headers_map
                .get("FamilyData")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            record_status: headers_map
                .get("RecordStatus")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            description: headers_map
                .get("Description")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            date_last_evaluated: headers_map
                .get("DateLastEvaluated")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
        };

        let key = (chr_norm, pos_parsed, ref_val.to_string(), alt_val.to_string());
        map.insert(key, rec);
    }

    println!("  -> Parsed {} lines from summary TSV.", line_count - 1);
    writeln!(
        log_file,
        "  -> Parsed {} lines from summary TSV.",
        line_count - 1
    )?;
    Ok(map)
}

/// Our final output record combining user input, ClinVar VCF, ClinVar TSV, and 1000G frequencies.
#[derive(Debug)]
struct FinalRecord {
    chr: String,
    pos: u32,
    ref_allele: String,
    alt_allele: String,
    clnsig: String,
    is_alt_pathogenic: bool,
    gene: Option<String>,
    allele_id: Option<i32>,
    genotype: String,
    review_stars: u8,
    af_esp: Option<f64>,
    af_exac: Option<f64>,
    af_tgp: Option<f64>,
    clndn: Option<String>,
    molecular_consequence: Option<String>,
    functional_consequence: Option<String>,
    mode_of_inheritance: Option<String>,
    preferred_values: Option<String>,
    citations: Option<String>,
    comments: Option<String>,
    family_data: Option<String>,
    record_status: Option<String>,
    description: Option<String>,
    date_last_evaluated: Option<String>,
    af_afr: Option<f64>,
    af_amr: Option<f64>,
    af_eas: Option<f64>,
    af_eur: Option<f64>,
    af_sas: Option<f64>,
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut log_file = OpenOptions::new()
        .create(true)
        .append(true)
        .open("pathogenic.log")?;

    println!("=== Pathogenic Variant Finder ===");
    writeln!(log_file, "=== Pathogenic Variant Finder ===")?;

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
        return Err("Genome build must be 'GRCh37' or 'GRCh38'".into());
    }
    if !input_path.exists() {
        return Err(format!("Input VCF not found: {}", input_path.display()).into());
    }

    // 1) CLINVAR VCF / TBI
    let (clinvar_url, clinvar_tbi_url) = match build.as_str() {
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
    fs::create_dir_all(&clinvar_dir)?;
    let clinvar_vcf_path = clinvar_dir.join(format!("clinvar_{}.vcf.gz", build));
    let clinvar_tbi_path = clinvar_dir.join(format!("clinvar_{}.vcf.gz.tbi", build));

    // 2) CLINVAR TSV for deeper annotation
    // We'll store it in the same directory
    let tsv_path = clinvar_dir.join("clinvar_summary_pathogenic.tsv.xz");
    let tsv_url = "https://raw.githubusercontent.com/SauersML/pathogenic/main/clinvar_summary_pathogenic.tsv.xz";

    // 3) 1000 GENOMES FREQUENCY FILE
    let onekg_url = match build.as_str() {
        "GRCH37" => "https://ftp.ensembl.org/pub/grch37/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz",
        "GRCH38" => "https://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz",
        _ => unreachable!(),
    };
    let onekg_file_path = clinvar_dir.join(format!("1000GENOMES-phase_3_{}.vcf.gz", build));

    println!("[STEP] Ensuring necessary files exist locally...");
    writeln!(log_file, "[STEP] Ensuring necessary files exist locally...")?;

    // Download CLINVAR VCF if not found or corrupt
    if !clinvar_vcf_path.exists() {
        println!("  -> Missing ClinVar VCF for {build}, downloading...");
        writeln!(
            log_file,
            "  -> Missing ClinVar VCF for {build}, downloading..."
        )?;
        download_file(clinvar_url, &clinvar_vcf_path, &mut log_file)?;
    } else {
        println!("  -> Found local ClinVar VCF: {}", clinvar_vcf_path.display());
        writeln!(log_file, "  -> Found local ClinVar VCF: {}", clinvar_vcf_path.display())?;
        // Check if corrupt
        let testf = File::open(&clinvar_vcf_path)?;
        let mut test_decoder = MultiGzDecoder::new(testf);
        let mut buffer = [0u8; 1024];
        if let Err(e) = test_decoder.read(&mut buffer) {
            println!("  -> ClinVar VCF is corrupt ({e}), re-downloading...");
            writeln!(
                log_file,
                "  -> ClinVar VCF is corrupt ({e}), re-downloading..."
            )?;
            fs::remove_file(&clinvar_vcf_path)?;
            download_file(clinvar_url, &clinvar_vcf_path, &mut log_file)?;
        }
    }

    // Download CLINVAR TBI if missing
    if !clinvar_tbi_path.exists() {
        println!("  -> Missing ClinVar TBI, downloading...");
        writeln!(log_file, "  -> Missing ClinVar TBI, downloading...")?;
        download_file(clinvar_tbi_url, &clinvar_tbi_path, &mut log_file)?;
    } else {
        println!("  -> Found local ClinVar TBI: {}", clinvar_tbi_path.display());
        writeln!(log_file, "  -> Found local ClinVar TBI: {}", clinvar_tbi_path.display())?;
    }

    // Download the 1000G frequency file if missing
    if !onekg_file_path.exists() {
        println!("  -> Missing 1000 Genomes frequency file for {build}, downloading...");
        writeln!(
            log_file,
            "  -> Missing 1000 Genomes frequency file for {build}, downloading..."
        )?;
        download_file(onekg_url, &onekg_file_path, &mut log_file)?;
    } else {
        println!("  -> Found local 1000 Genomes: {}", onekg_file_path.display());
        writeln!(log_file, "  -> Found local 1000 Genomes: {}", onekg_file_path.display())?;
        let testf = File::open(&onekg_file_path)?;
        let mut test_decoder = MultiGzDecoder::new(testf);
        let mut buffer = [0u8; 1024];
        if let Err(e) = test_decoder.read(&mut buffer) {
            println!("  -> 1000 Genomes freq file is corrupt ({e}), re-downloading...");
            writeln!(
                log_file,
                "  -> 1000 Genomes freq file is corrupt ({e}), re-downloading..."
            )?;
            fs::remove_file(&onekg_file_path)?;
            download_file(onekg_url, &onekg_file_path, &mut log_file)?;
        }
    }

    // Download ClinVar summary TSV if missing
    if !tsv_path.exists() {
        println!("  -> Missing ClinVar summary TSV, downloading...");
        writeln!(log_file, "  -> Missing ClinVar summary TSV, downloading...")?;
        download_file(tsv_url, &tsv_path, &mut log_file)?;
    } else {
        println!("  -> Found local ClinVar summary TSV: {}", tsv_path.display());
        writeln!(
            log_file,
            "  -> Found local ClinVar summary TSV: {}",
            tsv_path.display()
        )?;
    }

    // Parse the ClinVar VCF
    let (clinvar_map, clinvar_file_date) = parse_clinvar_vcf_gz(&clinvar_vcf_path, &mut log_file)?;
    println!("[LOG] ClinVar date: {}", clinvar_file_date);
    writeln!(log_file, "[LOG] ClinVar date: {}", clinvar_file_date)?;

    // Parse the ClinVar TSV for deeper annotation
    let tsv_map = parse_clinvar_tsv(&tsv_path, &build, &mut log_file)?;

    // Parse the user VCF input
    let input_variants = parse_input_vcf(&input_path, &mut log_file)?;

    println!("[STEP] Matching user variants with ClinVar and 1000G...");
    writeln!(
        log_file,
        "[STEP] Matching user variants with ClinVar and 1000G..."
    )?;

    // We'll store partial results first (from user input matched to ClinVar),
    // then combine with 1000G data, and finally combine with TSV info.
    #[derive(Debug)]
    struct TempRecord {
        chr: String,
        pos: u32,
        ref_allele: String,
        alt_allele: String,
        clnsig: String,
        is_alt_pathogenic: bool,
        gene: Option<String>,
        allele_id: Option<i32>,
        genotype: String,
        review_stars: u8,
        af_esp: Option<f64>,
        af_exac: Option<f64>,
        af_tgp: Option<f64>,
        clndn: Option<String>,
    }

    let pb = ProgressBar::new(input_variants.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=>-"),
    );

    let temp_results: Vec<TempRecord> = input_variants
        .par_iter()
        .flat_map_iter(|(_, iv)| {
            let mut local_found = Vec::new();
            for (alt_a, is_present) in &iv.alts {
                if !is_present {
                    continue;
                }
                let key = (iv.chr.clone(), iv.pos, iv.ref_allele.clone(), alt_a.clone());
                if let Some(cv) = clinvar_map.get(&key) {
                    if !cv.is_alt_pathogenic {
                        continue;
                    }
                    let genotype = iv.genotype.clone();
                    let review_stars = review_status_to_stars(cv.clnrevstat.as_deref());
                    local_found.push(TempRecord {
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
                        clndn: cv.clndn.clone(),
                    });
                }
            }
            pb.inc(1);
            local_found
        })
        .collect();

    pb.finish_with_message("ClinVar matching complete.");

    println!("  -> Matched {} variants in ClinVar", temp_results.len());
    // Collect keys of interest from temp_results for 1000 Genomes frequency lookup
    let keys_of_interest: HashSet<(String, u32, String, String)> = temp_results
        .iter()
        .map(|r| (r.chr.clone(), r.pos, r.ref_allele.clone(), r.alt_allele.clone()))
        .collect();

    // Stream 1000 Genomes VCF and extract frequencies for keys_of_interest
    println!("[STEP] Querying 1000 Genomes frequency VCF using CSI index: {}", onekg_file_path.display());
    writeln!(log_file, "[STEP] Querying 1000 Genomes frequency VCF using CSI index: {}", onekg_file_path.display())?;
    let onekg_index_path = clinvar_dir.join(format!("1000GENOMES-phase_3_{}_vcf.gz.csi", build));
    if !onekg_index_path.exists() {
        println!("  -> Missing 1000 Genomes CSI index for {}. Downloading...", build);
        writeln!(log_file, "  -> Missing 1000 Genomes CSI index for {}. Downloading...", build)?;
        download_file(
            if build == "GRCH37" {
                "https://ftp.ensembl.org/pub/grch37/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz.csi"
            } else {
                "https://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz.csi"
            },
            &onekg_index_path,
            &mut log_file
        )?;
    }

    let file = File::open(&onekg_file_path)?;
    let index_file = File::open(&onekg_index_path)?;
    let mut csi_reader = csi::io::Reader::new(index_file);
    let index = csi_reader.read_index()?;
    let mut reader = IndexedReader::new(file, index);
    let header = reader.read_header()?;

    let mut onekg_freqs: HashMap<(String, u32, String, String), (Option<f64>, Option<f64>, Option<f64>, Option<f64>, Option<f64>)> = HashMap::with_capacity(keys_of_interest.len());
    let mut unique_regions: HashSet<(String, u32)> = HashSet::new();
    for (chr, pos, _ref, _alt) in keys_of_interest.iter() {
        unique_regions.insert((chr.clone(), *pos));
    }

    for (chr, pos) in unique_regions {
        let region_str = format!("{}:{}-{}", chr, pos, pos);
        let region_obj: Region = region_str.parse()?;
        if let Ok(mut query) = reader.query(&header, &region_obj) {
            while let Some(record) = query.next() {
                let record = record?;
                let record_chr = record.reference_sequence_name().to_string();
                let record_pos = record
                    .variant_start()
                    .expect("Failed to get variant start")?;
                let record_ref = record.reference_bases().to_string();
                let record_alts: Vec<String> = record
                    .alternate_bases()
                    .iter()
                    .filter_map(Result::ok)
                    .map(|a| a.to_string())
                    .collect();
                let info = record.info();
                for (alt_idx, alt) in record_alts.iter().enumerate() {
                    let key = (record_chr.clone(), record_pos.get() as u32, record_ref.clone(), alt.clone());
                    if keys_of_interest.contains(&key) {
                        // Helper function to extract the frequency for a specific alternate allele from the INFO field
                        let get_freq = |key: &str, alt_idx: usize| -> Option<f64> {
                            match info.get(&header, key) {
                                Some(Ok(Some(Value::Array(array)))) => {
                                    match array {
                                        noodles::vcf::variant::record::info::field::value::Array::Float(values) => {
                                            let mut iter = values.iter();
                                            if let Some(Ok(Some(f))) = iter.nth(alt_idx) {
                                                Some(f as f64)
                                            } else {
                                                None
                                            }
                                        },
                                        _ => None, // Handle other Array variants (e.g., Integer, String)
                                    }
                                },
                                _ => None,
                            }
                        };
                        let afr = get_freq("AFR", alt_idx);
                        let amr = get_freq("AMR", alt_idx);
                        let eas = get_freq("EAS", alt_idx);
                        let eur = get_freq("EUR", alt_idx);
                        let sas = get_freq("SAS", alt_idx);
                        onekg_freqs.insert(key, (afr, amr, eas, eur, sas));
                    }
                }
            }
        }
    }
    println!("  -> Extracted frequencies for {} variants", onekg_freqs.len());
    writeln!(log_file, "  -> Extracted frequencies for {} variants", onekg_freqs.len())?;
    
    writeln!(
        log_file,
        "  -> Matched {} variants in ClinVar",
        temp_results.len()
    )?;

    // Next, we combine with 1000G frequency data
    let mut final_records: Vec<FinalRecord> = Vec::new();
    for r in &temp_results {
        let key = (r.chr.clone(), r.pos, r.ref_allele.clone(), r.alt_allele.clone());
        let (af_afr, af_amr, af_eas, af_eur, af_sas) = match onekg_freqs.get(&key) {
            None => (None, None, None, None, None),
            Some((afr, amr, eas, eur, sas)) => (*afr, *amr, *eas, *eur, *sas),
        };
        let annotation = tsv_map.get(&key);
        let final_rec = FinalRecord {
            chr: r.chr.clone(),
            pos: r.pos,
            ref_allele: r.ref_allele.clone(),
            alt_allele: r.alt_allele.clone(),
            clnsig: r.clnsig.clone(),
            is_alt_pathogenic: r.is_alt_pathogenic,
            gene: r.gene.clone(),
            allele_id: r.allele_id,
            genotype: r.genotype.clone(),
            review_stars: r.review_stars,
            af_esp: r.af_esp,
            af_exac: r.af_exac,
            af_tgp: r.af_tgp,
            clndn: r.clndn.clone(),
            molecular_consequence: annotation.and_then(|ann| ann.molecular_consequence.clone()),
            functional_consequence: annotation.and_then(|ann| ann.functional_consequence.clone()),
            mode_of_inheritance: annotation.and_then(|ann| ann.mode_of_inheritance.clone()),
            preferred_values: annotation.and_then(|ann| ann.preferred_values.clone()),
            citations: annotation.and_then(|ann| ann.citations.clone()),
            comments: annotation.and_then(|ann| ann.comments.clone()),
            family_data: annotation.and_then(|ann| ann.family_data.clone()),
            record_status: annotation.and_then(|ann| ann.record_status.clone()),
            description: annotation.and_then(|ann| ann.description.clone()),
            date_last_evaluated: annotation.and_then(|ann| ann.date_last_evaluated.clone()),
            af_afr,
            af_amr,
            af_eas,
            af_eur,
            af_sas,
        };
        final_records.push(final_rec);
    }

    // Sort final records by chromosome, then by position
    final_records.sort_by(|a, b| {
        let c = a.chr.cmp(&b.chr);
        if c == std::cmp::Ordering::Equal {
            a.pos.cmp(&b.pos)
        } else {
            c
        }
    });

    // Output CSV
    println!("[STEP] Writing final CSV to stdout...");
    writeln!(log_file, "[STEP] Writing final CSV to stdout...")?;
    let mut wtr = csv::Writer::from_writer(std::io::stdout());

    wtr.write_record(&[
        "Chromosome",
        "Position",
        "Reference Allele",
        "Alternate Allele",
        "Clinical Significance",
        "Is Alt Pathogenic",
        "Gene",
        "ClinVar Allele ID",
        "CLNDN",
        "Genotype",
        "Review Stars",
        "AF_ESP",
        "AF_EXAC",
        "AF_TGP",
        "Molecular Consequence",
        "Functional Consequence",
        "Mode of Inheritance",
        "Preferred Values",
        "Citations",
        "Comments",
        "Family Data",
        "Record Status",
        "Description",
        "Date Last Evaluated",
        "AF_AFR",
        "AF_AMR",
        "AF_EAS",
        "AF_EUR",
        "AF_SAS",
    ])?;

    for rec in &final_records {
        wtr.write_record(&[
            &rec.chr,
            &rec.pos.to_string(),
            &rec.ref_allele,
            &rec.alt_allele,
            &rec.clnsig,
            &rec.is_alt_pathogenic.to_string(),
            rec.gene.as_deref().unwrap_or(""),
            &rec.allele_id.map(|id| id.to_string()).unwrap_or_default(),
            rec.clndn.as_deref().unwrap_or(""),
            &rec.genotype,
            &rec.review_stars.to_string(),
            &rec.af_esp.map(|f| f.to_string()).unwrap_or_default(),
            &rec.af_exac.map(|f| f.to_string()).unwrap_or_default(),
            &rec.af_tgp.map(|f| f.to_string()).unwrap_or_default(),
            rec.molecular_consequence.as_deref().unwrap_or(""),
            rec.functional_consequence.as_deref().unwrap_or(""),
            rec.mode_of_inheritance.as_deref().unwrap_or(""),
            rec.preferred_values.as_deref().unwrap_or(""),
            rec.citations.as_deref().unwrap_or(""),
            rec.comments.as_deref().unwrap_or(""),
            rec.family_data.as_deref().unwrap_or(""),
            rec.record_status.as_deref().unwrap_or(""),
            rec.description.as_deref().unwrap_or(""),
            rec.date_last_evaluated.as_deref().unwrap_or(""),
            &rec.af_afr.map(|f| f.to_string()).unwrap_or_default(),
            &rec.af_amr.map(|f| f.to_string()).unwrap_or_default(),
            &rec.af_eas.map(|f| f.to_string()).unwrap_or_default(),
            &rec.af_eur.map(|f| f.to_string()).unwrap_or_default(),
            &rec.af_sas.map(|f| f.to_string()).unwrap_or_default(),
        ])?;
    }

    wtr.flush()?;
    println!(
        "Done. Wrote {} variants to CSV with ClinVar VCF data, ClinVar TSV data, and 1000G allele frequencies.",
        final_records.len()
    );
    writeln!(
        log_file,
        "Done. Wrote {} variants to CSV with ClinVar VCF data, ClinVar TSV data, and 1000G allele frequencies.",
        final_records.len()
    )?;

    Ok(())
}
