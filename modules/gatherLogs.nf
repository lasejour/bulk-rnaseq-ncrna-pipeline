#!/usr/bin/env nextflow

process GATHERLOGS {
        tag "Saving logs from all multiqc steps"

        memory '1 GB'
        publishDir "${params.output}/final_logs" , mode: 'copy'

        input:
        path LogOuts

        output:
        file LogOuts

        script:
        "echo printing"
}
