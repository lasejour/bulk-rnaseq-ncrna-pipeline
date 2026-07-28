#!/usr/bin/env nextflow

/*
Process for building a star index
Expects input channel to have:
1. A genome '.fa' file
2. An output directory to put the data
3. A '.gtf' file for genome annotation

Output:
1. The star index to  a results directory
*/

process STAR_INDEX {

        tag "building star index"

        publishDir "${params.output}/genomeIndices", mode: 'copy'

        input:
        tuple path(genome_file), path(genome_annot)
        path tx2gene_extraction

        output:
        path "starGenomeIndex", emit: star_index
        path "tx2gene_annotation.txt", emit: tx2gene_file

        script:
        """
        STAR --runThreadN $task.cpus --genomeDir starGenomeIndex \
             --runMode genomeGenerate --genomeFastaFiles $genome_file \
             --sjdbGTFfile $genome_annot

        Rscript $tx2gene_extraction $genome_annot
        """
}
