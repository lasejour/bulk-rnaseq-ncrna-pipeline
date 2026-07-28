#!/usr/bin/env nextflow

/*
Basic process for running mapping with Star
Expects input channel with:
- sample name, fastqs in tuple as first parameter
- path to genome index file as the second parameter
*/

// Mapping against genome with STAR
process STAR {

        tag "sample: $sample_name"

        publishDir "${params.output}/$sample_name/genomeMapping", mode: 'copy'

        input:
        tuple val(sample_name), path(trimmed_samples)
        path genome_index

        output:
        path "*"
        path "*Log.final.out", emit: star_logs
        tuple val(sample_name), path("*.Aligned.out.bam"), emit: bams

        script:
        """
        STAR --runThreadN $task.cpus --genomeDir $genome_index \
        --readFilesIn $trimmed_samples \
        --outFileNamePrefix ./$sample_name. \
        --limitBAMsortRAM 10000000000 \
        --outSJfilterDistToOtherSJmin 10 10 5 10 \
        --winAnchorMultimapNmax 100 \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outSAMtype BAM Unsorted \
        --outSJfilterCountUniqueMin 10 3 3 3 \
        --outSJfilterCountTotalMin 10 3 3 3 \
        --readFilesCommand zcat
        """
}
