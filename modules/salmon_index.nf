#!/usr/bin/env nextflow

/*
Process for building a salmon index
Expects input channel to have:
1. A transcriptome '.fa' file
2. A genome file
3. A corresponding genome annotation file
4. An output directory to put the data
5. The decoy transcriptome script

Output:
1. The salmon index to  a results directory
*/

process SALMON_INDEX {

        tag "building salmon index"

        publishDir "${params.output}/genomeIndices", mode: 'copy'

        input:
        path transcriptome_file
        path genome_file
        path genome_annot
        path decoy_script

        output:
        path "decoy_transcriptome"
        path "salmon_index", emit: sal_index

        script:
        """
        mashmap_cmd="\$(which mashmap)"
        bedtools_cmd="\$(which bedtools)"

        bash "${decoy_script}" -j ${task.cpus} -b "\${bedtools_cmd}" -m "\${mashmap_cmd}" \
              -a ${genome_annot} -g ${genome_file} -t ${transcriptome_file} \
              -o "decoy_transcriptome"

        awk '/^>/{sub(/ .*/,"")}1' decoy_transcriptome/gentrome.fa > decoy_transcriptome/new_gentrome.fa

        salmon index -p ${task.cpus} -t decoy_transcriptome/new_gentrome.fa -i "salmon_index" \
              --decoys decoy_transcriptome/decoys.txt -k 31

        """
}
