# WGBS pipeline
 This nextflow pipeline is used for processing whole genome bisulfite sequencing (WGBS) fastq data and has been used in the TOAST manuscript [ref]. 

# Processing steps
 This nextflow pipeline contains the following key processing steps: 
 1) Quality check using [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
 2) Split large fastq file into smaller data subsets.
 3) Trim adaptors and low quality bases using [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/). Subsequently, trim extra 18nt for Swiftbio libraries.
 4) Quality check (using [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) of trimmed reads from step 3.
 5) Map trimmed reads to GRCh38 reference genome using [Bismark algorithm](https://www.bioinformatics.babraham.ac.uk/projects/bismark/).
 6) Map trimmed reads to lambda DNA using [Bismark algorithm](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) to estimate bisulfite conversion rate.
 7) Merge processed data subsets (bam files) into a single file using [Samtools](http://www.htslib.org/).
 8) Deduplicate merged bam file using [Bismark algorithm](https://www.bioinformatics.babraham.ac.uk/projects/bismark/).
 9) Winnow out mapped reads with MAPQ score <5 using [Samtools](http://www.htslib.org/).
 10) Extract the methylation signal and calculate the methylation level using [Bismark algorithm](https://www.bioinformatics.babraham.ac.uk/projects/bismark/).
 11) Sort bam files (using [Samtools](http://www.htslib.org/)) for downstream analysis.
 12) Collect performance metrics using [Picard algorithm](https://broadinstitute.github.io/picard/).
 13) Determine the library insert size using [Picard algorithm](https://broadinstitute.github.io/picard/).

 Independent scripts available for the following downstream processing.
 1) bin/s1_bl_filtering.sh: To remove CpGs that reside in the [ENCODE blacklist](https://pubmed.ncbi.nlm.nih.gov/31249361/).
 2) bin/s2_extract_cov_gte5.sh: To extract CpGs with a sequencing coverage of at least 5x.
 3) bin/s3_SAS_filter.sh: To remove known (1000 Genome SAS) SNPs from the data.
<br/>

The canonical workflow for the WGBS pipeline is as below.

![wgbs](https://github.com/TOAST-LOLIPOP/WGBS-Pipeline/assets/143380802/ec896ded-3804-4904-bdd9-2f8ac39f92b3)

