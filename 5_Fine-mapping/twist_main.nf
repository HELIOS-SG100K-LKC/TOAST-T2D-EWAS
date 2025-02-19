
// IO Directories
prjDir="nf_twist"
rawFASTQ="RawFASTQ"
resultDir="twist_batchxx"

prjPath="/data/projects/"
resultPath="/data/projects/scratch/"

//HPC
userDir="/twist/"
twistSeqDir="nf_twist_pipeline"
toolsDir="/twist/"

params.outdir="${resultPath}/${userDir}/${twistSeqDir}/${prjDir}"

FASTQC_HOME="${prjPath}/${toolsDir}/opt/fastQC_v0.11.9"
TRIMGALORE_PATH="${prjPath}/${toolsDir}/opt/trimGalore_v0.6.4"
BWAMETH_PATH="${prjPath}/${toolsDir}/opt/bwaMeth_v0.2.2/bin"
GATK_PATH="${prjPath}/${toolsDir}/opt/gatk_v4.1.8.1"
METHYLDACKEL_PATH="${prjPath}/${toolsDir}/opt/methyldackel_v0.5.3"
PICARD_PATH="${prjPath}/${toolsDir}/opt/picard_v2.22.8"
SAMTOOLS_PATH="${prjPath}/${toolsDir}/opt/samtools_v1.10/samtools-1.10/bin"
SEQTK_PATH="${prjPath}/${toolsDir}/opt/seqtk_v1.3"
CUTADAPT_PATH="${prjPath}/${toolsDir}/opt/cutadapt_v3.2/bin/cutadapt"

GENOME_HG19_PATH="${prjPath}/${toolsDir}/ref_genome/hg19/"
BAIT_INTERVAL_PATH="${prjPath}/${toolsDir}/ref_genome/covered_target_bed/toast/"
TARGET_PATH="${prjPath}/${toolsDir}/ref_genome/covered_target_bed/toast/"

// the sample to be analyzed
read_pairs_ch = Channel.fromFilePairs("${baseDir}/$rawFASTQ/*_R{1,2}.fq.gz")
read_pairs_ch2 = Channel.fromFilePairs("${baseDir}/$rawFASTQ/*_R{1,2}.fq.gz")  


process raw_fastqc {

    tag "$sample_id"
    echo true
    
    publishDir "${params.outdir}/${resultDir}/${sample_id}/fastqc_raw", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) from read_pairs_ch

    output:
    tuple val(sample_id),file('*_fastqc.{zip,html}') into ch_fastqc_results_for_multiqc1
  
    script:
    """
    module load python/2.7.18

    ${FASTQC_HOME}/fastqc --noextract -t 2 ${reads[0]} ${reads[1]}

    """
}


process trim_galore {

    tag "$sample_id"
    echo true

    input:
     tuple val(sample_id), file(reads) from read_pairs_ch2

    output:
    tuple val(sample_id), file('*_R1.150bp_5prime.fq.gz'), file('*_R2.150bp_5prime.fq.gz') into ch_trimmed_reads_for_bwaMeth

    script:

    """             
    module rm python
    module load python/3.8.13

    ${TRIMGALORE_PATH}/trim_galore --cores 4 --gzip --path_to_cutadapt ${CUTADAPT_PATH} --hardtrim5 150 --2colour 20 --basename ${sample_id}_trimmed --paired ${reads[0]} ${reads[1]}
	    
    """
}


process bwaMeth {

    tag "$sample_id"
    echo true

    publishDir "${params.outdir}/${resultDir}/${sample_id}/bwaMeth", mode:'copy'

    input:
    tuple val(sample_id), file(reads_r1),file(reads_r2) from ch_trimmed_reads_for_bwaMeth

    output:
    tuple val(sample_id), file('*.sam') into ch_trimmed_reads_for_sambamba 
    tuple val(sample_id), file('*.sam') into ch_trimmed_reads_for_samtools_p1

    script:

    """

    module rm python
    module load python/3.8.13
    module load bwa
    module load sambamba

    ${BWAMETH_PATH}/bwameth.py --reference ${GENOME_HG19_PATH}/hg19.fa -t 12 --read-group '@RG\\tID:${sample_id}\\tPL:illumina\\tLB:${sample_id}\\tSM:${sample_id}' ${reads_r1} ${reads_r2} > ${sample_id}_v2.sam

    """
}


process sambamba {

    tag "$sample_id"
    echo true

    publishDir "${params.outdir}/${resultDir}/${sample_id}/bwaMeth", mode:'copy'

    input:
    tuple val(sample_id), file(sam_file) from ch_trimmed_reads_for_sambamba

    output:
    tuple val(sample_id), file('*_filtered_sorted.bam'), file('*_filtered_sorted.bam.bai') into ch_trimmed_reads_for_markDup, ch_trimmed_reads_for_samtools_p2

    script:

    """

    module load sambamba
    module load samtools

    sambamba view -h -t 4 -T ${GENOME_HG19_PATH}/hg19.fa --filter 'not secondary_alignment and not failed_quality_control and not supplementary and proper_pair and mapping_quality > 0' -f bam -l 0 -S ${sam_file} -o ${sample_id}_filtered.bam
    sambamba sort -t 4 -m 20GiB --tmpdir . -o /dev/stdout -l 0 ${sample_id}_filtered.bam | sambamba view -h -t 4 -o ${sample_id}_filtered_sorted.bam -T ${GENOME_HG19_PATH}/hg19.fa -f bam /dev/stdin
    samtools index -@ 4 ${sample_id}_filtered_sorted.bam > ${sample_id}_filtered_sorted.bam.bai

    """
}


process markDup {

    tag "$sample_id"
    echo true

    publishDir "${params.outdir}/${resultDir}/${sample_id}/markDup", mode:'copy'

    input:
    tuple val(sample_id), file(bam_file),file(bai_file) from ch_trimmed_reads_for_markDup

    output:
    tuple val(sample_id), file('*_filtered_sorted_markdup.bam'), file('*_filtered_sorted_markdup.bam.bai') into ch_trimmed_reads_for_methylDackel    
    tuple val(sample_id), file('*_filtered_sorted_markdup.bam'), file('*_filtered_sorted_markdup.bam.bai') into ch_trimmed_reads_for_picard
    tuple val(sample_id), file('*.txt') into ch_trimmed_reads_for_multiqc3

    script:
    
    """
    module load samtools    

    java -Xmx12g -Xms4g -jar ${PICARD_PATH}/picard.jar MarkDuplicates I=${bam_file} O=${sample_id}_filtered_sorted_markdup.bam R=${GENOME_HG19_PATH}/hg19.fa M=${sample_id}_picard_markdup_raw_metrics.txt CREATE_INDEX=false MAX_RECORDS_IN_RAM=1000 SORTING_COLLECTION_SIZE_RATIO=0.15 ASSUME_SORT_ORDER=coordinate OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 TMP_DIR=.
    samtools index -@ 4 ${sample_id}_filtered_sorted_markdup.bam

    """
}


process picard {

    tag "$sample_id"
    echo true


    publishDir "${params.outdir}/${resultDir}/${sample_id}/picard_after_markDup", mode:'copy'

    input:
    tuple val(sample_id), file(bam_file),file(bai_file) from ch_trimmed_reads_for_picard

    output:
    tuple val(sample_id), file('*') into ch_trimmed_reads_for_multiqc4


    script:

    """

    java -Xmx12g -Xms4g -jar ${PICARD_PATH}/picard.jar CollectHsMetrics I=${bam_file} O=${sample_id}_filtered_sorted_markdup_picard_collecthsmetrics_raw_metrics.txt R=${GENOME_HG19_PATH}/hg19.fa BAIT_INTERVALS=${BAIT_INTERVAL_PATH}/Probes_merged_ok_Methyl_JC_Sentinel_Ver5_MTE-99715360_hg19_230507132335.intervals TARGET_INTERVALS=${BAIT_INTERVAL_PATH}/Target_bases_covered_by_probes_Methyl_JC_Sentinel_Ver5_MTE-99715360_hg19_230507132335.intervals MINIMUM_MAPPING_QUALITY=20 COVERAGE_CAP=2500 PER_TARGET_COVERAGE=${sample_id}_filtered_sorted_markup_picard_collecthsmetrics_per_target_coverage.txt NEAR_DISTANCE=500 PER_BASE_COVERAGE=${sample_id}_filtered_sorted_markup_picard_collecthsmetrics_per_base_coverage.txt
    java -Xmx12g -Xms4g -jar ${PICARD_PATH}/picard.jar CollectMultipleMetrics I=${bam_file} O=${sample_id}_filtered_sorted_markup_picard_collectmultiplemetrics_raw.txt R=${GENOME_HG19_PATH}/hg19.fa PROGRAM=null PROGRAM=CollectGcBiasMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=CollectAlignmentSummaryMetrics

    """
}


process methylDackel {

    tag "$sample_id"
    echo true

    publishDir "${params.outdir}/${resultDir}/${sample_id}/methylDackel", mode:'copy'

    input:
    tuple val(sample_id), file(bam_file),file(bai_file) from ch_trimmed_reads_for_methylDackel

    output:
    tuple val(sample_id), file('*_CpG_processed.bedGraph') into ch_trimmed_reads_for_stats_cpg
    tuple val(sample_id), file('*_report.cytosine_report.txt') into ch_trimmed_reads_for_stats_non_cpg

    script:

    """
    
    ${METHYLDACKEL_PATH}/MethylDackel extract -@ 4 --keepDupes --minDepth 10 --maxVariantFrac 0.25 --OT 0,0,0,138 --OB 0,0,13,0 --mergeContext ${GENOME_HG19_PATH}/hg19.fa ${bam_file} -o ${sample_id}
    awk 'BEGIN {FS=OFS="\t"} NR == 1 {print \$0} NR > 1 {print \$1,\$2,\$3,((\$5/(\$5+\$6)*100)+0),\$5,\$6;}' OFMT="%.2f" ${sample_id}_CpG.bedGraph > ${sample_id}_CpG_processed.bedGraph
    ${METHYLDACKEL_PATH}/MethylDackel extract -@ 4 --keepDupes --minDepth 10 --maxVariantFrac 0.25 --OT 0,0,0,138 --OB 0,0,13,0 --cytosine_report --CHH --CHG ${GENOME_HG19_PATH}/hg19.fa ${bam_file} -o ${sample_id}_report

    """
}


