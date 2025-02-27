// src/main.rs

use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fmt;
use std::fs::{self, File};
use std::io::{copy, Read};
use std::path::{Path, PathBuf};

use clap::Parser;
use rayon::prelude::*;
use reqwest::blocking;
use rust_htslib::bcf::{header::HeaderView, Read as BcfRead, Reader};

/// Command-line arguments
#[derive(Parser)]
#[command(
    name = "pathogenic",
    about = "Identify pathogenic variants from a VCF using ClinVar data"
)]
struct Args {
    /// Genome build: "GRCh37" or "GRCh38"
    #[arg(short, long)]
    build: String,

    /// Path to the input VCF (uncompressed) file
    #[arg(short, long, value_name = "FILE")]
    input: PathBuf,
}

/// A small error type to unify download failures
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

/// Download a file from `url` to `out_path`.
fn download_file(url: &str, out_path: &Path) -> Result<(), DownloadError> {
    let response = blocking::get(url)?;
    let mut file = File::create(out_path)?;
    copy(&mut response.take(u64::MAX), &mut file)?;
    Ok(())
}

/// Helper to retrieve the name of the first contig in a BCF/VCF header
fn first_contig_name(header: &HeaderView) -> String {
    if header.contig_count() == 0 {
        return String::new();
    }
    // Attempt rid=0
    match header.rid2name(0) {
        Ok(name_bytes) => match std::str::from_utf8(name_bytes) {
            Ok(s) => s.to_string(),
            Err(_) => "".to_owned(),
        },
        Err(_) => "".to_owned(),
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    // Parse command-line arguments
    let args = Args::parse();
    let build = args.build.to_uppercase();
    let input_path = args.input.clone();

    // Validate genome build
    if build != "GRCH37" && build != "GRCH38" {
        eprintln!("Error: Genome build must be 'GRCh37' or 'GRCh38'.");
        std::process::exit(1);
    }
    // Validate input file
    if !input_path.exists() {
        eprintln!("Error: Input VCF not found: {}", input_path.display());
        std::process::exit(1);
    }

    // Decide ClinVar URLs
    let clinvar_dir = Path::new("clinvar_data");
    let clinvar_vcf_path = clinvar_dir.join(format!("clinvar_{}.vcf.gz", build));
    let clinvar_tbi_path = clinvar_vcf_path.with_extension("vcf.gz.tbi");
    let base_url = match build.as_str() {
        "GRCH37" => "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz",
        "GRCH38" => "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
        _ => unreachable!(),
    };
    let tbi_url = format!("{}.tbi", base_url);

    // Download ClinVar VCF + TBI if missing
    if !clinvar_vcf_path.exists() || !clinvar_tbi_path.exists() {
        fs::create_dir_all(&clinvar_dir)?;
        eprintln!("No local ClinVar data. Downloading for {build}...");

        let vcf_path_clone = clinvar_vcf_path.clone();
        let base_url_clone = base_url.to_string();
        let vcf_thread = std::thread::spawn(move || -> Result<(), DownloadError> {
            download_file(&base_url_clone, &vcf_path_clone)?;
            Ok(())
        });

        let tbi_path_clone = clinvar_tbi_path.clone();
        let tbi_url_clone = tbi_url.clone();
        let tbi_thread = std::thread::spawn(move || -> Result<(), DownloadError> {
            download_file(&tbi_url_clone, &tbi_path_clone)?;
            Ok(())
        });

        // Join threads
        match vcf_thread.join() {
            Ok(res) => {
                if let Err(e) = res {
                    eprintln!("Error downloading ClinVar VCF: {}", e);
                    std::process::exit(1);
                }
            }
            Err(e) => {
                eprintln!("Error: thread panicked downloading VCF: {:?}", e);
                std::process::exit(1);
            }
        }
        match tbi_thread.join() {
            Ok(res) => {
                if let Err(e) = res {
                    eprintln!("Error downloading ClinVar index: {}", e);
                    std::process::exit(1);
                }
            }
            Err(e) => {
                eprintln!("Error: thread panicked downloading TBI: {:?}", e);
                std::process::exit(1);
            }
        }
        eprintln!("ClinVar data downloaded to {}", clinvar_dir.display());
    }

    // Open ClinVar and input VCF
    let mut clinvar_reader = Reader::from_path(&clinvar_vcf_path)
        .map_err(|e| format!("Failed opening ClinVar VCF: {e}"))?;
    let mut input_reader = Reader::from_path(&input_path)
        .map_err(|e| format!("Failed opening input VCF: {e}"))?;

    // Determine "chr" usage
    let input_header = input_reader.header().clone();
    let clinvar_header = clinvar_reader.header().clone();
    let input_uses_chr = first_contig_name(&input_header).starts_with("chr");
    let clinvar_uses_chr = first_contig_name(&clinvar_header).starts_with("chr");

    /// Holds ClinVar data for a variant
    #[derive(Debug, Clone)]
    struct ClinVarInfo {
        clnsig: String,
        gene: Option<String>,
        allele_id: Option<i32>,
    }

    // (chr, pos, ref, alt) => ClinVarInfo
    let mut clinvar_variants: HashMap<(String, u32, String, String), ClinVarInfo> = HashMap::new();

    // ----------------- READ CLINVAR -----------------
    for rec_result in clinvar_reader.records() {
        let record = match rec_result {
            Ok(r) => r,
            Err(_) => continue,
        };
        let rid = match record.rid() {
            Some(r) => r,
            None => continue,
        };
        let chr_bytes = match clinvar_header.rid2name(rid) {
            Ok(b) => b,
            Err(_) => continue,
        };
        let mut chr = match std::str::from_utf8(chr_bytes) {
            Ok(s) => s.to_owned(),
            Err(_) => continue,
        };
        let pos = record.pos() + 1;

        // Normalize
        if input_uses_chr && !clinvar_uses_chr {
            if chr == "MT" {
                chr = "chrM".to_string();
            } else {
                chr = format!("chr{}", chr);
            }
        } else if !input_uses_chr && clinvar_uses_chr {
            if chr.eq_ignore_ascii_case("chrM") || chr.eq_ignore_ascii_case("chrMT") {
                chr = "MT".to_string();
            } else if let Some(stripped) = chr.strip_prefix("chr") {
                chr = stripped.to_string();
            }
        }

        // CLNSIG
        let clnsig_buf = match record.info(b"CLNSIG").string() {
            Ok(Some(buf)) => buf,
            _ => continue,
        };
        let clnsig_vals: Vec<String> = clnsig_buf
            .iter()
            .filter_map(|v| std::str::from_utf8(v).ok().map(|s| s.to_string()))
            .collect();
        if clnsig_vals.is_empty() {
            continue;
        }
        let clnsig_str = clnsig_vals.join("|");
        // Must contain "Pathogenic" or "Likely_pathogenic" but not "Conflicting..."
        if !clnsig_str.contains("Pathogenic")
            || clnsig_str.contains("Conflicting_interpretations_of_pathogenicity")
        {
            continue;
        }

        // GENEINFO => "SYMBOL:ID"
        let geneinfo = match record.info(b"GENEINFO").string() {
            Ok(Some(buf)) => {
                if let Some(raw) = buf.iter().next() {
                    if let Ok(s) = std::str::from_utf8(raw) {
                        Some(s.split(':').next().unwrap_or(s).to_string())
                    } else {
                        None
                    }
                } else {
                    None
                }
            }
            _ => None,
        };

        // ALLELEID => first integer if present
        let allele_id: Option<i32> = match record.info(b"ALLELEID").integer() {
            Ok(Some(bb)) => {
                // Typically there's 1 sub-array for an INFO with Number=1 or A.
                // We get the 0th array (if any) => &[i32], then the 0th i32
                bb.iter().next().copied()
            }
            _ => None,
        };

        let alleles = record.alleles();
        if alleles.is_empty() {
            continue;
        }
        let ref_allele = match std::str::from_utf8(alleles[0]) {
            Ok(a) => a.to_owned(),
            Err(_) => continue,
        };
        // Insert each ALT into map
        for alt_bytes in &alleles[1..] {
            if let Ok(alt_str) = std::str::from_utf8(alt_bytes) {
                let key = (chr.clone(), pos as u32, ref_allele.clone(), alt_str.to_owned());
                clinvar_variants.insert(
                    key,
                    ClinVarInfo {
                        clnsig: clnsig_str.clone(),
                        gene: geneinfo.clone(),
                        allele_id,
                    },
                );
            }
        }
    }

    // ----------------- PARSE INPUT VCF -----------------
    struct InputVariant {
        chr: String,
        pos: u32,
        ref_allele: String,
        alts: Vec<(String, bool)>, // alt allele + is present
    }
    let mut input_variants: Vec<InputVariant> = Vec::new();

    for rec_result in input_reader.records() {
        let record = match rec_result {
            Ok(r) => r,
            Err(_) => continue,
        };
        let rid = match record.rid() {
            Some(r) => r,
            None => continue,
        };
        let chr_bytes = match input_header.rid2name(rid) {
            Ok(b) => b,
            Err(_) => continue,
        };
        let mut chr = match std::str::from_utf8(chr_bytes) {
            Ok(s) => s.to_owned(),
            Err(_) => continue,
        };
        let pos = record.pos() + 1;
        // Normalize
        if input_uses_chr && !clinvar_uses_chr {
            if chr == "MT" {
                chr = "chrM".to_string();
            } else {
                chr = format!("chr{}", chr);
            }
        } else if !input_uses_chr && clinvar_uses_chr {
            if chr.eq_ignore_ascii_case("chrM") || chr.eq_ignore_ascii_case("chrMT") {
                chr = "MT".to_string();
            } else if let Some(stripped) = chr.strip_prefix("chr") {
                chr = stripped.to_string();
            }
        }

        let alleles = record.alleles();
        if alleles.is_empty() {
            continue;
        }
        let ref_allele = match std::str::from_utf8(alleles[0]) {
            Ok(a) => a.to_owned(),
            Err(_) => continue,
        };
        let alt_strings: Vec<String> = alleles[1..]
            .iter()
            .filter_map(|ab| std::str::from_utf8(ab).ok().map(|s| s.to_string()))
            .collect();

        // figure out which alt is present
        let mut present_flags = HashSet::new();
        match record.format(b"GT").integer() {
            Ok(gt_buf) => {
                // If there's at least one sample => check first sample
                if let Some(&gt_slice) = gt_buf.get(0) {
                    // gt_slice is &[i32]
                    for &allele_idx in gt_slice {
                        if allele_idx >= 1 {
                            present_flags.insert(allele_idx as usize);
                        }
                    }
                } else {
                    // no genotype for sample => assume all alt
                    for i in 1..alleles.len() {
                        present_flags.insert(i);
                    }
                }
            }
            Err(_) => {
                // parse error => assume all alt
                for i in 1..alleles.len() {
                    present_flags.insert(i);
                }
            }
        }

        let mut alt_pairs = Vec::new();
        for (i, alt_str) in alt_strings.into_iter().enumerate() {
            let alt_index = i + 1;
            let is_present = present_flags.contains(&alt_index);
            alt_pairs.push((alt_str, is_present));
        }

        input_variants.push(InputVariant {
            chr,
            pos: pos as u32,
            ref_allele,
            alts: alt_pairs,
        });
    }

    // ----------------- MATCH VS CLINVAR -----------------
    #[derive(Debug)]
    struct OutputRecord {
        chr: String,
        pos: u32,
        ref_allele: String,
        alt: String,
        clnsig: String,
        gene: Option<String>,
        allele_id: Option<i32>,
    }

    // Parallel lookups
    let mut results: Vec<OutputRecord> = input_variants
        .par_iter()
        .flat_map_iter(|var| {
            let mut found = Vec::new();
            for (alt, present) in &var.alts {
                if !present {
                    continue;
                }
                let key = (
                    var.chr.clone(),
                    var.pos,
                    var.ref_allele.clone(),
                    alt.clone(),
                );
                if let Some(info) = clinvar_variants.get(&key) {
                    found.push(OutputRecord {
                        chr: var.chr.clone(),
                        pos: var.pos,
                        ref_allele: var.ref_allele.clone(),
                        alt: alt.clone(),
                        clnsig: info.clnsig.clone(),
                        gene: info.gene.clone(),
                        allele_id: info.allele_id,
                    });
                }
            }
            found.into_iter()
        })
        .collect();

    // sort by chr then pos
    results.sort_by(|a, b| match a.chr.cmp(&b.chr) {
        std::cmp::Ordering::Equal => a.pos.cmp(&b.pos),
        o => o,
    });

    // ----------------- OUTPUT CSV -----------------
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
    for rec in results {
        wtr.write_record(&[
            rec.chr,
            rec.pos.to_string(),
            rec.ref_allele,
            rec.alt,
            rec.clnsig,
            rec.gene.unwrap_or_default(),
            rec.allele_id.map(|id| id.to_string()).unwrap_or_default(),
        ])?;
    }
    wtr.flush()?;
    Ok(())
}
