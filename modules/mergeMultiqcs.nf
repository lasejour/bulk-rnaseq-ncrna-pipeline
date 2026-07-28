#!/usr/bin/env nextflow

process MERGEMULTIQCS {
        tag "Merging all multiqc files into 1 file"

        memory '1 GB'
        publishDir "${params.output}/final_logs", mode: 'copy'

        input:
        path data
        path script

        output:
        file 'merged_results.tsv'

        script:
        """
        Rscript ${script} ${data}
        """
}
