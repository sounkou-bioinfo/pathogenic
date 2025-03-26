# Genetic Variant Finder

A high-performance Rust tool to identify pathogenic variants from VCF files using ClinVar and 1000 Genomes data.

## Overview

This tool analyzes VCF (Variant Call Format) files to identify potentially pathogenic genetic variants by cross-referencing them with the official ClinVar database and incorporating population frequency data from the 1000 Genomes Project. It supports both GRCh37 (hg19) and GRCh38 (hg38) genome builds.

## Features

- **Comprehensive Variant Detection**: Identifies variants classified as pathogenic, benign, or of uncertain significance in ClinVar
- **Population Frequency Integration**: Incorporates allele frequencies from 1000 Genomes (global and by population)
- **High Performance**:
  - Parallel processing using Rayon's work-stealing thread pool
  - Efficient file handling with memory-mapped files
  - Parallel chunk downloads for reference databases
- **Rich Annotation**: Provides comprehensive annotations including:
  - Clinical significance classifications
  - Gene information
  - Molecular and functional consequences
  - Mode of inheritance if available
  - Citations and additional evidence
- **Automatic Reference Data Management**: Downloads and maintains necessary reference databases
- **Advanced Reporting**: Generates detailed CSV reports and summary statistics files
- **Flexible Variant Selection**: Configurable inclusion of pathogenic, VUS, and benign variants

## Installation

```
# Clone the repository
git clone https://github.com/SauersML/pathogenic.git
cd pathogenic

# Build the project
cargo build --release

# The binary will be available at target/release/pathogenic
```

## Usage

```
# Basic usage (pathogenic variants only)
pathogenic --build GRCh38 --input your_variants.vcf

# Include Variants of Uncertain Significance (VUS)
pathogenic --build GRCh38 --input your_variants.vcf --include-vus

# Include benign variants
pathogenic --build GRCh38 --input your_variants.vcf --include-benign

# Include both VUS and benign variants
pathogenic --build GRCh38 --input your_variants.vcf --include-vus --include-benign
```

### Command-line Arguments

- `--build`, `-b`: Genome build, must be either "GRCh37" (hg19) or "GRCh38" (hg38)
- `--input`, `-i`: Path to the input VCF file (can be uncompressed or gzipped)
- `--include-vus`: Include variants of uncertain significance in the output
- `--include-benign`: Include benign variants in the output

## Output

The tool generates two output files in the `reports/` directory:

### 1. CSV Report File

The CSV report contains detailed information about all identified variants with the following columns:

- Chromosome, Position, Reference Allele, Alternate Allele
- Clinical Significance
- Is Alt Pathogenic
- Significance Category (pathogenic, benign, vus, conflicting)
- Gene
- ClinVar Allele ID
- Clinical Disease Name (CLNDN)
- Genotype
- Review Stars (0-4 star rating system)
- Allele frequencies:
  - ESP, ExAC, and 1000 Genomes Project (TGP)
  - Population-specific frequencies (AFR, AMR, EAS, EUR, SAS)
- Detailed annotations:
  - Molecular and Functional Consequences
  - Mode of Inheritance
  - Preferred Values
  - Citations
  - Comments and Family Data
  - Record Status
  - Description
  - Date Last Evaluated
 
Depending on available data, many of these fields may not be available.

### 2. Statistics Text File

A companion statistics file summarizing the analysis settings and results:

- Analysis settings (input file, genome build, variant types included)
- Total variants processed and reported
- Number of unique genes
- Counts of each variant classification type
- Distribution of allele frequencies across populations

### Output File Naming

Output files follow this naming convention:
```
[input_filename]_[analysis_type]_[timestamp].csv
[input_filename]_[analysis_type]_[timestamp]_stats.txt
```

Where `analysis_type` indicates which variant types were included in the report.

## How It Works

1. **Data Collection**: The tool automatically downloads necessary reference databases:
   - ClinVar VCF (GRCh37 or GRCh38)
   - ClinVar summary TSV for deeper annotation
   - 1000 Genomes frequency data

2. **Variant Processing**:
   - Parses the user's VCF input
   - Filters for variants present in the sample's genotype
   - Matches variants against ClinVar based on selected variant types
   - Integrates 1000 Genomes allele frequency data
   - Adds detailed annotations from ClinVar summary data

3. **Output Generation**:
   - Sorts variants by chromosome and position
   - Outputs comprehensive CSV with all annotations
   - Generates a statistics file summarizing the analysis

## Logging

The tool maintains a log file (`pathogenic.log`) that captures all processing steps and can be useful for troubleshooting.

## Documentation

For more detailed information, see the documentation in the `docs/` directory:

- [Implementation Details](docs/implementation_details.md)
- [Noodles Integration](docs/noodles_integration.md)
- [Parallel Processing](docs/parallel_processing.md)
- [Reporting Features](docs/reporting_features.md)
- [1000 Genome Frequency Extraction](docs/1000genome_frequency_extraction.md)

## Usage
Do not use this for clinical purposes.
- I might have written a bug in the code.
- Your input file might be wrong.
- Your input file may be right, but different somehow than the files I tested on.
- I did not test this extensively. In fact, there are no unit tests. (Though I did run an integration test via other tooling.)
- You might have done something wrong. For example, you entered the wrong genome build.
- The data in the database might be wrong about pathogenicity of some variant.
- You don't know how to interpret the results.
- You don't know how to interpret mode of inheritance.
- You don't know how to assess evidence for or confidence in pathogenicity.
- Your input file is imputed or has spurious genotype calls.
- Many other reasons not listed here.

What can you use this for?
- Finding potentially pathogenic variants to follow up on via other means.

## Disclaimer
This software is provided for educational purposes only. Results generated by the Genetic Variant Finder should not be used for medical diagnosis or to inform clinical decisions without proper validation by qualified healthcare professionals. The developers make no warranties, express or implied, regarding the accuracy, completeness, or reliability of the information provided by this tool. Users should be aware that variant classifications and pathogenicity assessments may change over time as new evidence emerges. Always consult with certified genetic counselors, clinical geneticists, or other appropriate healthcare providers for the interpretation of genetic variants in a medical context.
