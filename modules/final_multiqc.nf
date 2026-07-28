#!/usr/bin/env nextflow

process FINALMULTIQC {
        tag "Gathering logs..."

        publishDir "${params.output}/multiQC/finalMultiQC/", mode: 'copy'

        input:
        tuple val(GeneTypeStrand), path(file_names)
        path log_dir
        path config_file

        output:
        path '*'

        script:
        """
        multiqc $log_dir $file_names -o $GeneTypeStrand -c $config_file
        """

}
