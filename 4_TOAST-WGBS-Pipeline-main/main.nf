#!/usr/bin/env nextflow

// IO Directories
prjDir="wgbs/"
resultDir="wgbs_output"
rawFASTQ="RawFASTQ"

//on NSCC
userDir="darwin/"
wgbsDir="WGBS_2020/"
toolsDir="wgbs_tools/"

params.outdir="${userDir}/${wgbsDir}/${prjDir}"
FASTQC_HOME="${toolsDir}/opt/FastQC"
Q30="${toolsDir}/opt/q30"
TRIMGALORE_PATH="${toolsDir}/opt/TrimGalore-0.6.6"
BISMARK_PATH="${toolsDir}/opt/Bismark-0.22.3"
CUTADAPT_PATH="~/.local/bin/cutadapt"
BOWTIE_PATH="/app/bowtie2/2.29/"
GENOME_PATH="${toolsDir}/ref_genome/homo_sapiens/GRCh38"
GENOME_PRIMARY="${toolsDir}/ref_genome/homo_sapiens/GRCh38_primary"
GENOME_ALT="${toolsDir}/ref_genome/homo_sapiens/GRCh38_alt"
LAMBDA_PATH="${toolsDir}/ref_genome/lambda"
SAMTOOLS_PATH="${toolsDir}/opt/samtools/samtools-0.1.19/"
JAVA_PATH="${toolsDir}/opt/java/jdk-18/bin/"

PRESEQ_PATH="${toolsDir}/opt/preseq/preseq_v2.0"
QUALIMAP_PATH="${toolsDir}/opt/qualimap/qualimap_v2.2.1"
PICARD="${toolsDir}/opt/picard/picard.jar"
GENOME_FASTA="${toolsDir}/ref_genome/homo_sapiens/GRCh38/hg38.fa"


// Samples to be analyzed
read_pairs_ch = Channel.fromFilePairs("${baseDir}/$rawFASTQ/*_{1,2}.fastq.gz")


/******* STEP1 - FastQC_raw *********
***************************
*/
process fastqc_raw {

    tag "$sample_id"
    echo true
    
    publishDir "${params.outdir}/${resultDir}/${sample_id}/fastqc_raw", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) from read_pairs_ch

    output:
    tuple val(sample_id),file('*_fastqc.{zip,html}') into ch_fastqc_results_for_multiqc
    
    script:
    """
    module load python/2.7.12

    ${FASTQC_HOME}/fastqc -t 2 ${reads[0]} ${reads[1]}

    """
}


/******* STEP2 - fastq_split *********
***************************
*/

process fastq_split {

    tag "$sample_id"
    echo true
    
    input:
    tuple val(sample_id), file(reads_r1), file(reads_r2) from read_pairs_ch2

    output:
    tuple val(sample_id), file('*_1.*'), file('*_2.*') into ch_split_fastq

    script:
    """
    zcat $reads_r1 | split -l 80000000 -d - ${sample_id}_1. & 
    zcat $reads_r2 | split -l 80000000 -d - ${sample_id}_2. &
    wait

    """

}

// regroup the tuple by using two keys
ch_split_fastq.transpose()
        .map {input -> tuple(input[0], input[1].toString().tokenize('.').get(1), input[1], input[2])}
        .groupTuple(by: [0,1])
        .into { ch_fastq_for_trim1; ch_fastq_for_trim2; ch_fastq_for_trim3 }


/******* STEP3 - Trim Galore! *********
***************************
*/
process trim_galore {

    tag "$sample_id"
    echo true
 
    publishDir "${params.outdir}/${resultDir}/${sample_id}/trimmed_reads", pattern: "*trimming_report.txt", mode:'copy'
 
   input:
    tuple val(sample_id), val(sub_ID), file(split_reads_r1),file(split_reads_r2) from ch_fastq_for_trim3

    output:
    tuple val(sample_id),val(sub_ID), file('*_val_1.fq'), file('*val_2.fq') into ch_trimmed_reads_for_alignment,ch_trimmed_reads_fastQC
    
    tuple val(sample_id), file("*trimming_report.txt") into ch_trim_galore_results_for_multiqc

    script:

    """             
    module rm python
    module load python/3.5.1

    $TRIMGALORE_PATH/trim_galore --cores 4 --path_to_cutadapt ${CUTADAPT_PATH} --clip_R2 18 --three_prime_clip_R1 18 --phred33 --paired ${split_reads_r1} ${split_reads_r2}
    rm $split_reads_r1 $split_reads_r2
	    
    """
}

ch_trimmed_reads_for_alignment
	.groupTuple(by:[0,1], sort:true)
	.into { ch_trimmed_reads_for_align_genome; ch_trimmed_reads_for_align_lambda}

ch_trimmed_reads_fastQC
	.map {input -> tuple(input[0], input[2], input[3])}
	.groupTuple(by:0)
	.into {ch_trimmed_reads_fastQC2; ch_q30_aftertrim_reads1; ch_q30_aftertrim_reads2 }


/******* STEP4.1 - FastQC for trimmed reads *********
***************************
*/
process fastqc_aftertrim {

    tag "$sample_id"
    echo true
    
    publishDir "${params.outdir}/${resultDir}/${sample_id}/fastqc_aftertrim", pattern:"*_fastqc.{zip,html}", mode: 'copy'

    input:
    tuple val(sample_id), file(trimmed_reads_r1),file(trimmed_reads_r2) from ch_trimmed_reads_fastQC2

    output:
    tuple val(sample_id), file('*_fastqc.{zip,html}') into ch_fastqc_aftertrim_results_for_multiqc
    
    script:

    """
    module load python/2.7.12
    
    cat ${trimmed_reads_r1} > ${sample_id}_aftertrim_R1.fq &
    cat ${trimmed_reads_r2} > ${sample_id}_aftertrim_R2.fq &
    wait

    
    $FASTQC_HOME/fastqc -t 2 ${sample_id}_aftertrim_R1.fq ${sample_id}_aftertrim_R2.fq 
    rm ${sample_id}_aftertrim_R1.fq ${sample_id}_aftertrim_R2.fq

    """
}


/******* STEP4.2 - mapping to lambda with Bismark *********
**************************
*/
process bismark_align_lambda {

    tag "$sample_id"
    echo true
    
    publishDir "${params.outdir}/${resultDir}/${sample_id}/bam_files_lambda", mode:"copy"

    input:
    tuple val(sample_id), val(sub_ID), file(trimmed_reads_r1),file(trimmed_reads_r2) from ch_trimmed_reads_for_align_lambda

    output:
    tuple val(sample_id),file("*_bismark_bt2_PE_report.txt") into ch_bismark_lambda_for_multiqc

    script:
    """    
    module load bowtie2/2.29 
    module load samtools/1.3

    $BISMARK_PATH/bismark -multicore 2 --bowtie2 -p 2 --bam --score_min L,0,-0.2 --path_to_bowtie $BOWTIE_PATH $LAMBDA_PATH -1 ${trimmed_reads_r1} -2 ${trimmed_reads_r2} 
    rm ${sample_id}*_bismark_bt2_pe.bam
	   
    """
}


/******* STEP4.3 - mapping with Bismark *********
***************************
*/
process bismark_align {

    tag "$sample_id"
    echo true

    publishDir "${params.outdir}/${resultDir}/${sample_id}/bam_files", pattern:"*_bismark_bt2_PE_report.txt", mode: "copy"

    input:
    tuple val(sample_id), val(sub_ID), file(trimmed_reads_r1),file(trimmed_reads_r2) from ch_trimmed_reads_for_align_genome

    output:
    tuple val(sample_id), val(sub_ID), file('*_val_1_bismark_bt2_pe.bam'), file('*_val_1.fq_unmapped_reads_1_bismark_bt2_pe.bam') into ch_bismark_hg19_bam
    tuple val(sample_id), file("*_bismark_bt2_PE_report.txt") into ch_bismark_hg19_for_multiqc

    script:
    """   
    module load bowtie2/2.29 
    module load samtools/1.3

    ## map the overall reads to the hg38 primary assembly
    $BISMARK_PATH/bismark -multicore 4 --bowtie2 -p 3 --bam --un --ambiguous --score_min L,0,-0.2 --path_to_bowtie $BOWTIE_PATH $GENOME_PRIMARY -1 ${trimmed_reads_r1} -2 ${trimmed_reads_r2}
    
    
    ## map the unmapped reads to hg38 alternate contigs
    $BISMARK_PATH/bismark -multicore 4 --bowtie2 -p 3 --bam --score_min L,0,-0.2 --path_to_bowtie $BOWTIE_PATH $GENOME_ALT -1 ${trimmed_reads_r1}_unmapped_reads_1.fq.gz -2 ${trimmed_reads_r2}_unmapped_reads_2.fq.gz
    
    
    ## remove the unmapped reads and ambiguous reads
    rm ${sample_id}*_unmapped_reads_2.fq.gz ${sample_id}*_unmapped_reads_1.fq.gz & 
    rm ${sample_id}*_ambiguous_reads_2.fq.gz ${sample_id}*_ambiguous_reads_1.fq.gz &
    wait

    """
}


/******* STEP5 - merge bam files *********
***************************
*/
ch_bismark_hg19_bam
	.map{ input -> tuple(input[0], input[2], input[3]) }
	.groupTuple(by:0)
	.set {ch_bismark_hg19_bam2}

process merge_bam {

    tag "$sample_id"
    echo true
    
    publishDir "${params.outdir}/${resultDir}/${sample_id}/unsortedButMerged_ForBismark_file", mode: 'move'

    input:
    tuple val(sample_id), file(bam_file1), file(bam_file2) from ch_bismark_hg19_bam2

    output:    
    tuple val(sample_id), file("*_unsorted_merged.bam") into ch_bismark_hg19_merged_bam   

    script:
    """    
    module load samtools/1.3

    ## merge primary reads and alternate contig reads	
    samtools merge -nf -@ 6 ${sample_id}_unsorted_merged.bam $bam_file1 $bam_file2

    """
}


/******* STEP6 - deduplicate bam files *********
***************************
*/
process bam_dedup {
    tag "$sample_id"
    echo true
   
    publishDir "${params.outdir}/${resultDir}/${sample_id}/dedup_bams", pattern: "*deduplication_report.txt", mode:"copy"

    input:
    tuple val(sample_id), path(bamFile) from ch_bismark_hg19_merged_bam

    output:
    tuple val(sample_id), file('*.deduplicated.bam') into ch_deduplicated_bam_for_mapq
    tuple val(sample_id), file("*deduplication_report.txt") 

    script:
    """    
    module load samtools    

    $BISMARK_PATH/deduplicate_bismark -p --bam ${bamFile[0]} -o ${sample_id}
    

    """
}


/******* STEP7 - MAPQ5 filtering *********
***************************
*/
process mapq {
    tag "$sample_id"
    echo true
   
    input:
    tuple val(sample_id), path(bamFile) from ch_deduplicated_bam_for_mapq

    output:
    tuple val(sample_id), file('*_dedup_mapq5.bam') into ch_dedup_mapq5_bam_for_methylextract, ch_deduplicated_bam_for_sort
    

    script:
    """    
    module load samtools    
    
    samtools view -h -q 5 -@ 4 -o ${sample_id}_dedup_mapq5.bam ${sample_id}.deduplicated.bam

    """
}


/******* STEP8 - methylation extraction *********
***************************
*/
process methyl_extract {
    tag "$sample_id"
    echo true
   
    publishDir "${params.outdir}/${resultDir}/${sample_id}/methylation_extraction", mode:"copy"

    input:
    tuple val(sample_id), path(merged_dedup_mapq5_bam) from ch_dedup_mapq5_bam_for_methylextract

    output:
    file "*"

    script:
    """    
    $BISMARK_PATH/bismark_methylation_extractor -p --multicore 5 --gzip --no_overlap --comprehensive --merge_non_CpG --cutoff 1 --buffer_size 15G --zero_based --cytosine_report --genome_folder $GENOME_PATH $merged_dedup_mapq5_bam

    """

}


/******* STEP9 - sort bam *********
***************************
*/
process sort_bam {

    tag "$sample_id"
    echo true
    
    publishDir "${params.outdir}/${resultDir}/${sample_id}/sorted_bam", mode:"copy"
    
    input:
    tuple val(sample_id), file(merged_deduplicated_bam) from ch_deduplicated_bam_for_sort

    output:
    tuple val(sample_id), file('*sorted_deduplicated.bam') into ch_bismark_hg19_deduplicated_bam1, ch_bismark_hg19_deduplicated_bam2
   
    script:
    """
    samtools sort -@ 3 -m 3G $merged_deduplicated_bam ${sample_id}_sorted_deduplicated

    """

}


/******* STEP10 - coverage analysis *********
 ***************************
*/
process coverage {
    tag "$sample_id"
    echo true
    
    publishDir "${params.outdir}/${resultDir}/${sample_id}/coverage_files_picard",  mode:"copy"

    input:
    tuple val(sample_id), path(sorted_deduplicated_bam) from ch_bismark_hg19_deduplicated_bam1

    output:
    
    tuple val(sample_id),file("*_coverage_picard.txt")

    script:
    """   
    
    ${JAVA_PATH}/java -Xmx3G -jar $PICARD CollectWgsMetrics REFERENCE_SEQUENCE=$GENOME_FASTA \
 		MINIMUM_MAPPING_QUALITY=0 \
 		INPUT=$sorted_deduplicated_bam \
 		OUTPUT=${sample_id}_coverage_picard.txt

    """
}


/******* STEP11 - insert size analysis *********
* ***************************
*/
process insert_size {

    tag "$sample_id"
    echo true
    
    publishDir "${params.outdir}/${resultDir}/${sample_id}/insert_size_files", mode:"copy"

    input:
    tuple val(sample_id), path(sorted_deduplicated_bam) from ch_bismark_hg19_deduplicated_bam2

    output:
   
    tuple val(sample_id),file("*_insert_size_metrics.txt") 

    script:
    """    
    ${JAVA_PATH}/java -Xmx3G -jar $PICARD CollectInsertSizeMetrics \
        I=$sorted_deduplicated_bam \
        O=${sample_id}_insert_size_metrics.txt \
        H=${sample_id}_insert_size_histogram.pdf \
        INCLUDE_DUPLICATES=false \
        ASSUME_SORTED=true
	
    """

}

