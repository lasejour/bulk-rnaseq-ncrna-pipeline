#!/usr/bin/env nextflow


process ALTERNATIVESTAR {

        tag "sample: $sample_name"
        maxForks 12

        // Similar as previous steps, just added a few parameters to adjust the filename to make logs identifible
        publishDir "${params.output}/$sample_name/$tag", mode: 'copy'

        input:
        tuple val(sample_name), path(trimmed_samples)
        path genome_index
        val tag

        output:
        path "*Log.final.out", emit: star_logs

        script:
        """
        STAR --runThreadN $task.cpus --genomeDir $genome_index \
        --readFilesIn $trimmed_samples \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix ./$sample_name. \
        --limitBAMsortRAM 10000000000 \
        --seedPerWindowNmax 18 \
        --readFilesCommand zcat
        """
}
