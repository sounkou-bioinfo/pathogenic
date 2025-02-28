// src/main.rs

use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fmt;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, copy, Read};
use std::path::{Path, PathBuf};
use std::thread;

use clap::Parser;
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use reqwest::blocking;

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

/// A struct to hold ClinVar data for a variant
#[derive(Debug, Clone)]
struct ClinVarInfo {
    clnsig: String,
    gene: Option<String>,
    allele_id: Option<i32>,
}

/// Read the ClinVar VCF (which is .vcf.gz) line by line, parse out only the
/// Pathogenic/Likely_pathogenic variants, store them in a map keyed by (chr, pos, ref, alt).
fn parse_clinvar_vcf_gz(
    vcf_gz_path: &Path,
    input_uses_chr: bool,
) -> Result<HashMap<(String, u32, String, String), ClinVarInfo>, Box<dyn Error>> {
    let file = File::open(vcf_gz_path)?;
    let decoder = MultiGzDecoder::new(file);
    let reader = BufReader::new(decoder);

    let mut clinvar_variants = HashMap::new();
    let mut clinvar_uses_chr = false; // We'll guess after we see first contig in header

    let mut in_header = true;
    for line_res in reader.lines() {
        let line = match line_res {
            Ok(l) => l,
            Err(_) => continue,
        };
        // Skip blank lines
        if line.trim().is_empty() {
            continue;
        }
        if in_header {
            if line.starts_with("##") {
                // still in meta lines
                continue;
            } else if line.starts_with("#CHROM") {
                // header line
                in_header = false;
                continue;
            } else {
                // unexpected line, possibly #CHROM but missing the label
                continue;
            }
        }

        // Now we parse a data line
        // Format: CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  [FORMAT] [samples...]
        let mut fields = line.split('\t');
        let chrom_str = match fields.next() {
            Some(c) => c.to_string(),
            None => continue,
        };
        let pos_str = match fields.next() {
            Some(p) => p,
            None => continue,
        };
        // skip ID
        let _ = fields.next();
        let ref_str = match fields.next() {
            Some(r) => r,
            None => continue,
        };
        let alt_str = match fields.next() {
            Some(a) => a,
            None => continue,
        };
        // skip QUAL, FILTER
        let _ = fields.next();
        let _ = fields.next();
        // read INFO
        let info_str = match fields.next() {
            Some(i) => i,
            None => continue,
        };

        // Try to parse pos
        let pos_num = match pos_str.parse::<u32>() {
            Ok(pn) => pn,
            Err(_) => continue,
        };

        // Let's see if we guessed whether ClinVar has 'chr' or not
        // We'll do it after the first data line
        if clinvar_uses_chr == false {
            if chrom_str.starts_with("chr") {
                clinvar_uses_chr = true;
            }
        }

        // We'll unify naming with the input's style:
        // If input has chr and ClinVar doesn't, we add "chr" if it's not "MT" (or fix for "MT")
        // If input doesn't have chr and ClinVar does, we remove "chr"
        // We'll do that in a small function:
        let mut chr_fixed = chrom_str.clone();
        if input_uses_chr && !clinvar_uses_chr {
            if chr_fixed == "MT" {
                chr_fixed = "chrM".to_string();
            } else {
                chr_fixed = format!("chr{}", chr_fixed);
            }
        } else if !input_uses_chr && clinvar_uses_chr {
            if chr_fixed.eq_ignore_ascii_case("chrM") || chr_fixed.eq_ignore_ascii_case("chrMT") {
                chr_fixed = "MT".to_string();
            } else if let Some(stripped) = chr_fixed.strip_prefix("chr") {
                chr_fixed = stripped.to_string();
            }
        }

        // Parse multiple ALTs if any
        let alt_list: Vec<&str> = alt_str.split(',').collect();

        // Parse the INFO field into a map
        let mut info_map = parse_info_field(info_str);

        // Check CLNSIG
        let clnsig_opt = info_map.remove("CLNSIG"); // e.g. "Pathogenic|otherstuff"
        if clnsig_opt.is_none() {
            continue;
        }
        let clnsig_str = clnsig_opt.unwrap();
        // If not "Pathogenic" or "Likely_pathogenic", skip
        // Also skip if "Conflicting_interpretations_of_pathogenicity"
        if !clnsig_str.contains("Pathogenic") {
            continue;
        }
        if clnsig_str.contains("Conflicting_interpretations_of_pathogenicity") {
            continue;
        }

        // geneinfo
        let gene_opt = info_map.remove("GENEINFO").and_then(|val| {
            // e.g. "BRCA2:675"
            let part = val.split(':').next().unwrap_or(&val);
            Some(part.to_string())
        });

        // allele_id
        let allele_id_opt = info_map.remove("ALLELEID").and_then(|val| {
            // we might parse it as i32
            if let Ok(num) = val.parse::<i32>() {
                Some(num)
            } else {
                None
            }
        });

        // For each ALT, create an entry
        let ref_allele = ref_str.to_string();
        for alt_item in alt_list {
            let alt_allele = alt_item.to_string();

            let key = (chr_fixed.clone(), pos_num, ref_allele.clone(), alt_allele);
            let info_struct = ClinVarInfo {
                clnsig: clnsig_str.clone(),
                gene: gene_opt.clone(),
                allele_id: allele_id_opt,
            };
            clinvar_variants.insert(key, info_struct);
        }
    }

    Ok(clinvar_variants)
}

/// Parse the user's local, uncompressed VCF line by line
/// Return a vector of input variants, each with the alt alleles and whether they
/// are present in genotype (we do a simplistic approach).
fn parse_input_vcf(
    input_vcf_path: &Path,
    need_chr: bool,
) -> Result<Vec<InputVariant>, Box<dyn Error>> {
    let file = File::open(input_vcf_path)?;
    let reader = BufReader::new(file);

    let mut variants = Vec::new();
    let mut in_header = true;
    let mut input_uses_chr = false;

    // We'll determine if input has "chr" from first data line's CHROM
    // But we also have "need_chr" and "has_chr_clinvar" to unify naming
    for line_res in reader.lines() {
        let line = match line_res {
            Ok(l) => l,
            Err(_) => continue,
        };
        if line.trim().is_empty() {
            continue;
        }
        if in_header {
            if line.starts_with("##") {
                continue;
            } else if line.starts_with("#CHROM") {
                in_header = false;
                continue;
            } else {
                // skip other header lines
                continue;
            }
        }

        // parse data line
        let mut cols = line.split('\t');
        let chrom_str = match cols.next() {
            Some(c) => c.to_string(),
            None => continue,
        };
        let pos_str = match cols.next() {
            Some(c) => c,
            None => continue,
        };
        // skip ID
        let _ = cols.next();
        let ref_str = match cols.next() {
            Some(r) => r,
            None => continue,
        };
        let alt_str = match cols.next() {
            Some(a) => a,
            None => continue,
        };
        // skip QUAL, FILTER
        let _ = cols.next();
        let _ = cols.next();
        // read INFO
        let _ = cols.next().unwrap_or("");

        // read FORMAT + genotype if any
        // or skip if not present
        let format_and_more: Vec<&str> = cols.collect();

        // detect if input has "chr"
        if !input_uses_chr {
            if chrom_str.starts_with("chr") {
                input_uses_chr = true;
            }
        }

        // unify naming
        let mut chr_fixed = chrom_str.clone();
        if need_chr && !input_uses_chr {
            if chr_fixed == "MT" {
                chr_fixed = "chrM".to_string();
            } else {
                chr_fixed = format!("chr{}", chr_fixed);
            }
        } else if !need_chr && input_uses_chr {
            if chr_fixed.eq_ignore_ascii_case("chrM") || chr_fixed.eq_ignore_ascii_case("chrMT") {
                chr_fixed = "MT".to_string();
            } else if let Some(stripped) = chr_fixed.strip_prefix("chr") {
                chr_fixed = stripped.to_string();
            }
        }

        let pos_num = match pos_str.parse::<u32>() {
            Ok(pn) => pn,
            Err(_) => continue,
        };
        let ref_allele = ref_str.to_string();

        let alt_list: Vec<String> = alt_str.split(',').map(|s| s.to_string()).collect();

        // simplistic genotype parse: check if there's a "GT" column
        // in the format. If not, we'll assume all alt are present.
        let mut present_flags = HashSet::new();
        // format_and_more[0] might be FORMAT, subsequent are sample fields
        if !format_and_more.is_empty() {
            let format_cols = format_and_more[0].split(':').collect::<Vec<&str>>();
            // find which index is GT if any
            let gt_index_opt = format_cols.iter().position(|f| *f == "GT");

            if gt_index_opt.is_some() && format_and_more.len() > 1 {
                // we have at least one sample
                let sample_fields = format_and_more[1]; // first sample
                let sample_items = sample_fields.split(':').collect::<Vec<&str>>();
                let gt_index = gt_index_opt.unwrap();
                if gt_index < sample_items.len() {
                    let gt_val = sample_items[gt_index];
                    // e.g. "0/1" or "1|1"
                    // parse out numeric allele indices
                    let separators = ['/', '|'];
                    let alleles_gt = gt_val.split(&separators[..]).collect::<Vec<&str>>();
                    for al_str in alleles_gt {
                        if let Ok(idx) = al_str.parse::<usize>() {
                            if idx >= 1 {
                                // alt allele
                                present_flags.insert(idx);
                            }
                        }
                    }
                } else {
                    // no GT in that sample => assume all alt
                    for i in 1..=alt_list.len() {
                        present_flags.insert(i);
                    }
                }
            } else {
                // no GT field => assume all alt
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

        let mut alts_with_pres = Vec::new();
        for (i, alt_a) in alt_list.into_iter().enumerate() {
            let alt_idx = i + 1;
            let is_present = present_flags.contains(&alt_idx);
            alts_with_pres.push((alt_a, is_present));
        }

        variants.push(InputVariant {
            chr: chr_fixed,
            pos: pos_num,
            ref_allele,
            alts: alts_with_pres,
        });
    }

    Ok(variants)
}

/// A struct to represent one input VCF variant record and which alt alleles are present
#[derive(Debug)]
struct InputVariant {
    chr: String,
    pos: u32,
    ref_allele: String,
    alts: Vec<(String, bool)>,
}

/// A helper function to parse INFO field (semicolon-delimited key=value pairs)
fn parse_info_field(info: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for item in info.split(';') {
        if item.is_empty() {
            continue;
        }
        // might be KEY=VAL or just KEY
        let mut eqsplit = item.splitn(2, '=');
        let key = eqsplit.next().unwrap_or("").to_string();
        let val = eqsplit.next().unwrap_or("").to_string();
        // store
        map.insert(key, val);
    }
    map
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
        let vcf_thread = thread::spawn(move || -> Result<(), DownloadError> {
            download_file(&base_url_clone, &vcf_path_clone)?;
            Ok(())
        });

        let tbi_path_clone = clinvar_tbi_path.clone();
        let tbi_url_clone = tbi_url.clone();
        let tbi_thread = thread::spawn(move || -> Result<(), DownloadError> {
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

    // Parse the ClinVar VCF
    // We'll parse it as .vcf.gz with line-based approach
    // Then see if we need chr prefix or not after we see the data
    // We'll guess that the user wants chr prefix if input has it
    let user_temp_file = File::open(&input_path)?;
    let user_buf = BufReader::new(user_temp_file);
    let mut user_has_chr = false;
    // read lines until we see the first data line
    // then see if CHROM has 'chr'
    let mut in_header = true;
    for line_res in user_buf.lines() {
        let line = match line_res {
            Ok(l) => l,
            Err(_) => break,
        };
        if line.trim().is_empty() {
            continue;
        }
        if in_header {
            if line.starts_with("##") {
                continue;
            } else if line.starts_with("#CHROM") {
                in_header = false;
                continue;
            } else {
                continue;
            }
        } else {
            // data line
            let chrom_part = line.split('\t').next().unwrap_or("");
            if chrom_part.starts_with("chr") {
                user_has_chr = true;
            }
            break;
        }
    }

    // parse ClinVar with knowledge of user_has_chr
    let clinvar_map = parse_clinvar_vcf_gz(&clinvar_vcf_path, user_has_chr)?;

    // parse user input
    let input_variants = parse_input_vcf(&input_path, user_has_chr)?;

    // Now match them
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

    let mut results: Vec<OutputRecord> = input_variants
        .par_iter()
        .flat_map_iter(|var| {
            let mut found = Vec::new();
            for (alt, is_present) in &var.alts {
                if !is_present {
                    continue;
                }
                let key = (var.chr.clone(), var.pos, var.ref_allele.clone(), alt.clone());
                if let Some(info) = clinvar_map.get(&key) {
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

    // sort
    results.sort_by(|a, b| {
        let c = a.chr.cmp(&b.chr);
        if c == std::cmp::Ordering::Equal {
            a.pos.cmp(&b.pos)
        } else {
            c
        }
    });

    // output CSV
    let mut wtr = csv::Writer::from_writer(std::io::stdout());
    wtr.write_record(&["Chromosome","Position","Ref","Alt","ClinicalSignificance","Gene","ClinVarAlleleID"])?;
    for rec in results {
        wtr.write_record(&[
            rec.chr,
            rec.pos.to_string(),
            rec.ref_allele,
            rec.alt,
            rec.clnsig,
            rec.gene.clone().unwrap_or_default(),
            rec.allele_id.map(|id| id.to_string()).unwrap_or_default()
        ])?;
    }
    wtr.flush()?;
    Ok(())
}
