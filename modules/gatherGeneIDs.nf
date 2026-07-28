#!/usr/bin/env nextflow

process GATHERGENEIDS {
        tag "collecting salmon counts"

        memory '1 GB'
        publishDir "${params.output}", mode: 'copy'

        input:
        path trans_gene
        path tx2gene_salmon
        path salmon_quant_files

        output:
        file 'gene_expressions_values.tsv'

        script:
        """
        Rscript ${tx2gene_salmon} "." ${trans_gene}
        """
}
