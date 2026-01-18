To improve your repository's accessibility and usability for the research community, I recommend a structured `README.md` that reflects the multi-stage nature of your pipelines (R, Nextflow, and Shell).

Based on the repository structure, here is a template you can adapt.

---

# TOAST-T2D-EWAS

This repository contains the analytical pipelines and code for the study: **"Nuclear regulatory disturbances precede and predict the development of Type-2 diabetes in Asian populations"** (Jain, Ng, Tay, et al.).

Reference: [https://doi.org/10.1101/2025.02.14.25322264](https://doi.org/10.1101/2025.02.14.25322264)

## Table of Contents

* [System Requirements](https://www.google.com/search?q=%23system-requirements)
* [Installation Guide](https://www.google.com/search?q=%23installation-guide)
* [Repository Structure](https://www.google.com/search?q=%23repository-structure)
* [How to Run](https://www.google.com/search?q=%23how-to-run)
* [Data Preparation](https://www.google.com/search?q=%23data-preparation)
* [License](https://www.google.com/search?q=%23license)

---

## System Requirements

### Hardware

* **Memory:** Minimum 16GB RAM for small-scale analysis; 64GB+ recommended for WGBS and EWAS pipelines.
* **Storage:** High-performance storage for large genomic (VCF/BAM) and methylation (BigWig/BedGraph) files.

### Software & Dependencies

* **Operating System:** Linux (CentOS/Ubuntu recommended) or macOS.
* **Core Languages:** - R (>= 4.0.0)
* Python (>= 3.8)


* **Workflow Management:** [Nextflow](https://www.google.com/search?q=https://www.nextflow.io/) (for WGBS and RNA-seq pipelines).
* **Containerization (Optional but Recommended):** Docker or Singularity.

### R Dependencies

Ensure the following packages are installed:

```r
install.packages(c("tidyverse", "data.table", "BiocManager"))
BiocManager::install(c("limma", "minfi", "TOAST", "DMRcate", "coloc", "susieR"))

```

---

## Installation Guide

1. **Clone the repository:**
```bash
git clone https://github.com/HELIOS-SG100K-LKC/TOAST-T2D-EWAS.git
cd TOAST-T2D-EWAS

```


2. **Set up the environment:**
We recommend using Conda to manage environments:
```bash
conda env create -f environment.yml # If you provide one
conda activate toast-t2d

```


3. **Nextflow setup:**
Ensure Nextflow is executable in your path:
```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow

```



---

## Repository Structure

The code is organized by the chronological order of analysis:

* `1_EWAS-Pipeline`: Scripts for Epigenome-Wide Association Studies.
* `2_Functional-Annotation-and-Enrichment`: Post-EWAS annotation and pathway analysis.
* `3_Genomic-Analysis`: Integration with genetic variants (mQTLs).
* `4_TOAST-WGBS-Pipeline-main`: Nextflow pipeline for Whole Genome Bisulfite Sequencing.
* `5_Fine-mapping`: SuSiE-based fine-mapping of associations.
* `6_Methylation-Risk-Scores`: Generation and validation of MRS.
* `7 - RNA-seq pipeline`: Transcriptomic integration.
* `Data-Blobs`: Helper scripts and reference datasets.

---

## How to Run

### 1. Running the EWAS Pipeline (R)

Navigate to the EWAS directory and run the main analysis script:

```bash
Rscript 1_EWAS-Pipeline/run_ewas_analysis.R --input metadata.csv --output results/

```

### 2. Running the WGBS Pipeline (Nextflow)

This pipeline handles raw data processing. Edit the `nextflow.config` file to match your cluster environment (e.g., SLURM, SGE).

```bash
nextflow run 4_TOAST-WGBS-Pipeline-main/main.nf -profile <docker/singularity/conda> --input samplesheet.csv

```

### 3. Fine-mapping & Risk Scores

Follow the specific RMarkdown or scripts within directories `5_` and `6_` for specific downstream statistical modeling.

---

## Data Preparation

* **Methylation Data:** Expects IDAT files or pre-processed Beta-value matrices.
* **Metadata:** Ensure your phenotype file contains required columns: `Sample_ID`, `T2D_Status`, `Age`, `Sex`, and `Cell_Type_Proportions`.

---

## License

This project is licensed under the **GPL-3.0 License** - see the [LICENSE](https://www.google.com/search?q=LICENSE) file for details.

### Citation

If you use this code, please cite:

> Jain, A., Ng, T. K., Tay, W. T., et al. "Nuclear regulatory disturbances precede and predict the development of Type-2 diabetes in Asian populations." *medRxiv* (2025).

---

### Key Improvements Made:

1. **Logical Flow:** Organized the analysis from raw data processing (WGBS) to high-level analysis (MRS).
2. **Explicit Tooling:** Mentioned **Nextflow** and **R Bioconductor** dependencies, which are critical for bioinformaticians.
3. **Actionable Commands:** Added shell blocks for cloning and running scripts so users can "copy-paste."
4. **Hardware Warnings:** Added a memory warning since EWAS and WGBS data are notoriously memory-intensive.
