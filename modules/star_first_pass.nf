#!/usr/bin/env nextflow

/*
STAR splicing process - First Pass
Discovers splice junctions

Expects input channel with:
- sample name, fastqs in tuple as first parameter
- path to genome index file as the second parameter
It outputs:
- Splice junction table
*/

process STAR_FIRST_PASS {

    tag "sample: $sample_name"
    conda 'bioconda::star=2.7.11b'


    input:
    tuple val(sample_name), path(read_files)
    path genome_index

    output:
    tuple val(sample_name), path("*SJ.out.tab"), emit: SJ_logs

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
        --outFileNamePrefix ./$sample_name. \\
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
        --outSAMtype None \\
        --outSJfilterCountUniqueMin 10 3 3 3 \\
        --outSJfilterCountTotalMin 10 3 3 3 \\
        --readFilesCommand zcat
    """
}
