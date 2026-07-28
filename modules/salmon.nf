#!/usr/bin/env nextflow

/*
Basic process for running Quantification with salmon
Expects input channel with:
- sample name, fastqs in tuple as first parameter
- path to salmon index as the second parameter
- file for a perl script to determine strandedness of fastq files
- output directory path as third parameter
*/

// Mapping against transcriptome with salmon
process SALMON {

        tag "sample: $sample_name"
        maxForks 12

        publishDir "${params.output}/$sample_name/salmonResults", mode: 'copy'

        input:
        tuple val(sample_name), path(read_files)
        path salmon_index
        path stranded_script

        output:
        path "*", emit: sal_logs
        stdout emit: strand_value
        path "${sample_name}/${sample_name}_quant.sf", emit: quant_file

        script:
        if(params.prep == "single-end") {
              """
              salmon quant --threads $task.cpus --validateMappings --libType=A -i $salmon_index -r $read_files -o $sample_name --seqBias --gcBias
              perl $stranded_script "${sample_name}/logs/salmon_quant.log"
              mv "${sample_name}/quant.sf" "${sample_name}/${sample_name}_quant.sf"
              """
        } else {
              """
              salmon quant --threads $task.cpus --validateMappings --libType=A -i $salmon_index -1 ${read_files[0]} -2 ${read_files[1]} -o $sample_name --seqBias --gcBias
              perl $stranded_script "${sample_name}/logs/salmon_quant.log"
              mv "${sample_name}/quant.sf" "${sample_name}/${sample_name}_quant.sf"
              """
        }

}
