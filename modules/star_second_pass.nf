#!/usr/bin/env nextflow

/*
STAR splicing process - second Pass
Aligns reads using cohort level merged junction

Expects input channel with:
- sample name, fastqs in tuple as first parameter
- path to genome index file as the second parameter
- path merged splice junction table

It outputs:
- Summary of mapping
- Splice junction per sample
- BAM file
*/

process STAR_SECOND_PASS {

    tag "sample: $sample_name"
    conda 'bioconda::star=2.7.11b'

    publishDir "${params.output}/$sample_name/2nd_genomeMapping/", mode: 'copy'

    input:
    tuple val(sample_name), path(read_files)
    path genome_index
    path merged_sj

    output:
    tuple val(sample_name), path("*.Aligned.out.bam"), emit: bams
    path "*.Log.final.out", emit: star_logs
    path "*.SJ.out.tab", emit: sj_pass
    path "*.ReadsPerGene.out.tab"

    script:
    def reads_arg = read_files instanceof List && read_files.size() > 1
        ? "${read_files[0]} ${read_files[1]}"
        : "${read_files[0]}"

    """
    STAR \\
        --runMode alignReads \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${genome_index} \\
        --readFilesIn ${reads_arg} \\
        --sjdbFileChrStartEnd ${merged_sj} \\
        --outFileNamePrefix ./$sample_name. \\
        --quantMode GeneCounts \\
        --limitBAMsortRAM 10000000000 \\
        --outSJfilterDistToOtherSJmin 10 10 5 10 \\
        --winAnchorMultimapNmax 100 \\
        --outFilterMultimapNmax 20 \\
        --alignSJoverhangMin 8 \\
        --alignSJDBoverhangMin 1 \\
        --outFilterMismatchNmax 999 \\
        --outFilterMismatchNoverReadLmax 0.1 \\
        --alignIntronMin 20 \\
        --alignIntronMax 1000000 \\
        --alignMatesGapMax 1000000 \\
        --outSAMtype BAM Unsorted \\
        --outSJfilterCountUniqueMin 10 3 3 3 \\
        --outSJfilterCountTotalMin 10 3 3 3 \\
        --readFilesCommand zcat
    """
}
