#!/usr/bin/env nextflow

process FEATURECOUNTSMULTIQC {
        tag "featureCounts MultiQC"

        publishDir "${params.output}/multiQC/featureCounts", mode: 'copy'

        input:
        tuple val(stranded), path(file_names)
        path config_file

        output:
        path "exonStranded/multiqc_data/*multiqc_featureCounts.txt", emit: exonStranded_out, optional: true
        path "exonUnstranded/multiqc_data/*multiqc_featureCounts.txt", emit: exonUnstranded_out, optional: true
        path "geneStranded/multiqc_data/*multiqc_featureCounts.txt", emit: geneStranded_out, optional: true
        path "geneUnstranded/multiqc_data/*multiqc_featureCounts.txt", emit: geneUnstranded_out, optional: true

        script:
        """
        multiqc $file_names -o ${stranded} -c $config_file
        mv "${stranded}/multiqc_data/multiqc_featurecounts.txt" "${stranded}/multiqc_data/${stranded}_multiqc_featureCounts.txt"
        """

}
