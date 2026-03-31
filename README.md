# gut-microbiome-dynamics

**High-resolution tracking of gut microbiome community dynamics during antibiotic perturbation in mice, using barcoded *E. coli* lineage tracing and 16S amplicon sequencing.**

---

## Overview

This project reanalyses a controlled mouse gut microbiome experiment in which barcoded *Escherichia coli* was introduced into the gut of antibiotic-treated mice. The barcoding system tags individual *E. coli* lineages with unique genetic identifiers, allowing us to track hundreds of thousands of bacterial lineages simultaneously — far beyond what conventional 16S sequencing can resolve.

The combination of lineage-level barcode data with family-level 16S community profiles creates a dual-resolution view of gut ecology: how individual *E. coli* lineages compete, drift, and go extinct, while the broader microbial community restructures around them under antibiotic pressure.

An unexpected finding — the bloom of *Paenibacillus macerans* to >50% relative abundance in colonised mice — led to a deep genomics investigation using Oxford Nanopore long-read sequencing of isolates, revealing that the organism's spectinomycin resistance is intrinsic (a ribosomal protein deletion), not acquired via horizontal gene transfer as originally suspected.

## Experimental Design

| Group | Mice | Treatment | Timepoints | Sequencing |
|-------|------|-----------|------------|------------|
| Colonised cohort 1 | m1–m4 | Barcoded *E. coli* + spectinomycin | 19 (3h – 16d) | 16S amplicon + barcode amplicon |
| Colonised cohort 2 | m5–m8 | Barcoded *E. coli* + spectinomycin | 18 (3h – 14d) | 16S amplicon + barcode amplicon |
| Antibiotic-only controls | c_m1–c_m4 | Spectinomycin only (no *E. coli*) | 10 (1d – 9d) | 16S amplicon only |

**Nanopore isolates** (from cohort 2 faecal samples):

| Isolate | Mouse | Coverage | Contigs | Species | Status |
|---------|-------|----------|---------|---------|--------|
| M5 | m5 | 103x | 1 (circular) | *P. macerans* | Primary reference |
| M6 | m6 | 101x | 1 (circular) | *P. macerans* | Clonal with M5 |
| M8 | m8 | 84x | 1 (circular) | *P. macerans* | Clonal with M5 |
| M7 | m7 | 9x | 7 | Mixed | Excluded (contaminated) |

## Analysis Pipeline

The analysis is organised into five sequential R modules plus a genomics investigation in Python/shell:

### Module 1 — 16S Community Composition (`analysis/01_16S/`)
Processes 16S amplicon data into family- and genus-level relative abundance tables. Generates stacked area plots and log-scale line plots showing community restructuring over time. Rarefaction thresholds set per cohort based on sequencing depth.

### Module 2 — Barcode Lineage Dynamics (`analysis/02_barcode/`)
Imports barcoded *E. coli* frequency data (>500,000 unique barcodes across 8 mice), filters technical noise, and tracks lineage dynamics over time. Includes:
- Lineage frequency trajectories (area and line plots)
- Hill diversity indices (richness, Shannon, Simpson)
- Effective population size (Ne) estimation via Wright-Fisher drift model
- Neutral drift testing
- Cross-cohort frequency spectrum comparison

### Module 3 — Barcode Clustering (`analysis/03_clustering/`)
Hierarchical clustering of barcode time series using Shape-Based Distance (SBD) to identify co-behaving lineage groups. Applies z-normalisation and average linkage. Clusters represent lineages with similar colonisation dynamics (e.g., early colonisers vs late bloomers vs transients).

### Module 4 — Co-clustering (`analysis/04_coclustering/`)
Jointly clusters barcode lineage groups with 16S family abundance time series using SBD. This reveals which bacterial families co-vary with which *E. coli* lineage clusters — connecting strain-level dynamics to community-level ecology. Bootstrap validation via pvclust (1,000 replicates).

### Module 5 — Dynamic Community Modelling (`analysis/05_DCM/`)
Estimates time-varying Jacobian interaction matrices between all tracked species/lineage groups using expanding-window covariance. Extracts:
- Eigenvalue evolution (community stability over time)
- KPCA trajectories (state-space visualisation)
- Changepoint detection (when community stability shifts)
- Interaction strength outliers (which species pairs drive instability)
- Cross-mouse common outlier pairs (reproducible ecological interactions)

Key DCM finding: *Paenibacillaceae* is a structural outlier — C3/C4 barcode clusters interact with *Paenibacillaceae* in 7–8 of 12 mice.

### Genomics Investigation (`analysis/paenibacillus_resistance/`)
Nanopore long-read sequencing + comparative genomics to determine the mechanism of spectinomycin resistance in *P. macerans* isolates. Tested 10 hypotheses; resolved all. See `results/genomics/reports/paenibacillus_resistance_report.md` for the full narrative.

## Key Findings

**Community dynamics:**
- Spectinomycin treatment causes rapid community restructuring within 24–48 hours
- *Paenibacillaceae* blooms to 50–67% relative abundance in colonised mice but remains <1.4% in controls
- *E. coli* lineage diversity collapses over time consistent with Wright-Fisher drift + selection
- Barcode cluster–16S family co-clustering reveals reproducible ecological associations across mice

**Spectinomycin resistance mechanism:**
- *P. macerans* isolates carry an intrinsic 3-amino-acid deletion in ribosomal protein S5 (rpsE): **Δ(A21-K22-V23)**
- The deleted lysine (K22) is conserved across *E. coli*, *B. subtilis*, and the *P. macerans* type strain — its absence disrupts the spectinomycin binding site
- No acquired resistance genes (aadA, ANT(9), aminoglycoside-modifying enzymes) found anywhere in the genome
- Two integrative conjugative elements (108 kb and 25 kb) were acquired but carry no resistance cargo
- Resistance is **intrinsic**, not horizontally transferred — the original HGT hypothesis is not supported

## Repository Structure

```
gut-microbiome-dynamics/
├── analysis/
│   ├── 01_16S/              # 16S community composition
│   ├── 02_barcode/           # Barcode lineage dynamics + popgen
│   ├── 03_clustering/        # SBD hierarchical clustering
│   ├── 04_coclustering/      # Joint barcode–16S clustering
│   ├── 05_DCM/               # Jacobian stability analysis
│   └── paenibacillus_resistance/  # Genomics scripts (Python)
├── data/                     # Raw data (not tracked — see Data section)
│   ├── raw/                  # 16S tables, barcode frequencies
│   ├── references/           # Reference genomes
│   └── sequence_data/        # Nanopore FASTQ + assemblies
├── metadata/
│   └── 0_config.R            # Global config: samples, paths, palettes, themes
├── results/
│   ├── figures/              # All plots (16S, barcodes, clustering, DCM, genomics)
│   ├── tables/               # All tabular outputs
│   └── genomics/             # Genomic analysis outputs + reports
├── renv.lock                 # R package lockfile
└── CLAUDE.md                 # Development notes
```

## Getting Started

### Prerequisites

**R analysis (modules 1–5):**
- R >= 4.3
- Package management via [renv](https://rstudio.github.io/renv/) — run `renv::restore()` to install all dependencies from the lockfile
- Key packages: tidyverse, patchwork, here, TSclust, pvclust, changepoint, kernlab

**Genomics analysis:**
- Python >= 3.10 with matplotlib, seaborn, pandas, numpy, biopython
- WSL (Windows Subsystem for Linux) with:
  - minimap2 2.30+
  - samtools 1.23+
  - MAFFT v7.526+
  - Bakta v1.12+
  - Panaroo 1.5+
  - fastANI 1.34+
  - MUMmer (nucmer) 4.0+

### Running the pipeline

Each module sources its dependencies. From the project root:

```r
# Run modules in order
source("analysis/01_16S/04_plot_community.R")   # sources 00–03 internally
source("analysis/02_barcode/03_dynamics.R")       # sources 00–02 internally
source("analysis/03_clustering/02_select.R")      # sources 00–01 internally
source("analysis/04_coclustering/01_coclustering.R")
source("analysis/05_DCM/05_outliers.R")           # sources 01–04 internally
```

All paths use `here::here()` — no `setwd()` required.

## Data

Raw data files are **not tracked in git** (`.gitignore` excludes `data/`). The `data/` directory must be populated locally before running the pipeline. This includes:
- 16S amplicon abundance tables (`data/raw/16S/`)
- Barcode frequency matrices (`data/raw/barcodes/`)
- Nanopore FASTQ and assembly files (`data/sequence_data/`)
- Reference genomes (`data/references/`)

## Author

**Melis Gencel**

## License

Please contact the author before reuse.
