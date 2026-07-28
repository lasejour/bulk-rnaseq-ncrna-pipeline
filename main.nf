#!/usr/bin/env nextflow

include { FASTQC } from './modules/fastqc.nf'
include { FASTQC as cleanedQC } from './modules/fastqc.nf'

include { CUTADAPT } from './modules/cutadapt.nf'

include { MULTIQC } from './modules/multiqc.nf'
include { MULTIQC as cleanedMULTIQC } from './modules/multiqc.nf'
include { MULTIQC as salmonMULTIQC } from './modules/multiqc.nf'
include { MULTIQC as starMULTIQC } from './modules/multiqc.nf'
include { MULTIQC as riboMULTIQC } from './modules/multiqc.nf'
include { MULTIQC as spacerMULTIQC } from './modules/multiqc.nf'
include { FINALMULTIQC as finalMULTIQC_STRANDED } from './modules/final_multiqc.nf'
include { FINALMULTIQC as finalMULTIQC_UNSTRANDED } from './modules/final_multiqc.nf'

include { FEATURECOUNTSMULTIQC } from './modules/featureCounts_multiqc.nf'

include { SALMON } from './modules/salmon.nf'
include { STAR } from './modules/star.nf'
include { ALTERNATIVESTAR as RIBOSTAR } from './modules/alternative_star.nf'
include { ALTERNATIVESTAR as SPACERSTAR } from './modules/alternative_star.nf'

include { STAR_INDEX } from './modules/star_index.nf'
include { SALMON_INDEX } from './modules/salmon_index.nf'

include { READCOUNTS } from './modules/featureCounts.nf'

include { GATHERGENEIDS } from './modules/gatherGeneIDs.nf'
include { GATHERLOGS } from './modules/gatherLogs.nf'
include { MERGEMULTIQCS } from './modules/mergeMultiqcs.nf'


include { STAR_FIRST_PASS } from './modules/star_first_pass.nf'
include { MERGE_SJ } from './modules/merge_sj.nf'
include { STAR_SECOND_PASS } from './modules/star_second_pass.nf'

def helpMessage() {
    log.info """
    =========================================
     RNA-Seq Pipeline ncRNA core v2.0
    =========================================
    Usage:
        nextflow run RNA-Seq-pipeline-ncRNA-core.nf [options]

    Required arguments:
    * Either
        --fastq_dir               Path to input FASTQ files (e.g. 'data/')
        --pattern                 Regex to find FASTQ files (e.g. '*{R1,R2}*.fastq.gz')

    * Or
        --sample_file             Path to sample list with sample names and read locations (i.e refer to README.md for format)


    Optional arguments:
        --output                  Path to results directory (default: 'results')
        --prep                    Value of Read Orientation (either 'single-end' or 'paired-end')
        --genome_index            Path to directory of STAR genome index (default: defined in the params.yml)
        --genome_id_file          Path to file with ensembl transcript and gene ids (default: defined in the params.yml)
        --genome_annot_file       Path to genome gtf file (default: defined in the params.yml)

        --salmon_index            Path to directory of Salmon transcriptome index (default: defined in the params.yml)

        --rRNA_index              Path to directory ribosomal genome (default: defined in the params.yml)
        --rRNA_Spacers_index      Path to directory of ribosomal spacers genome (default: defined in the params.yml)

        --splice_junctions        Will splice junctions be accounted for? (either true or false)
        --min_unique_reads        Minimum number of reads for uniquely mapping reads across all sample (default: 3)

        * If you are constructing a STAR index:
        --genome_fasta_file       Path to primary assembly fasta file
        --genome_annot_file       Path to genome gtf file (default: defined in the params.yml)
        --genome_id_file          Path to file with ensembl transcript and gene ids (default: defined in the params.yml)

        * If you are constructing a salmon index:
        --transcriptome_file      Path to transcriptome fasta file

        * If you'd like to run on conda envirnoment
        -profile                 (e.g local)

        * Defining your own parameters
        -params-file             Path to parameters file (default: params.yml)


    Example:
        nextflow run RNA-Seq-pipeline-ncRNA-core.nf \\
            --fastq_dir 'data/' \\
            --pattern '*_{1,2}.fastq.gz'
            --output my_results
            -profile local
            -params-file params.yml
    =========================================
    """.stripIndent()

}

params.help = ""

if (params.help) {
    helpMessage()
    exit 0
}

// Genome requirements
if (!params.genome_index) {
    log.info "Because no genome index was given, STAR will construct an index.\n\n"
    if(!params.genome_fasta_file || !params.genome_annot_file) {
        exit 1, """Please provide the necessary file for STAR to construct an index:
        --genome_fasta_file
        --genome_annot_file"""
    }
}

// Transcriptome requirements
if (!params.salmon_index) {
    log.info "Because no salmon index was given, Salmon will construct an index.\n\n"
    if( !params.genome_fasta_file || !params.genome_annot_file || !params.transcriptome_file) {
        exit 1, """Please provide the necessary files for salmon to construct an index.
        They are:
        --transcriptome_file
        --genome_fasta_file
        --genome_annot_file"""
    }
}

workflow {
      // Check inputs

      // Prep
      if(params.prep == "single-end") {

        if (params.sample_file) {
            Channel
                  .fromPath(params.sample_file, checkIfExists: true)
                  .splitCsv(header: true, strip: true)
                  .map { row ->
                        def sample = row.SampleName
                        def r1     = file(row.Path, checkIfExists: true)

                        if (!sample) error "Missing 'SampleName' field in row: ${row}"
                        if (!r1)     error "Path not found for sample: ${sample}"
                        return [sample, r1]
                  }
                  .set { raw_reads_ch }

        } else if (!params.sample_file && params.fastq_dir && params.pattern ) {
            Channel
                  .fromPath("${params.fastq_dir}/${params.pattern}", checkIfExists: true)
                  .ifEmpty { exit 1, "You provided an empty directory" }
                  .map { file -> tuple(file.simpleName, file) }
                  .set { raw_reads_ch }

        } else {
            error """
              No input provided. Please specify either:
                --input      Path to a CSV samplesheet OR
                --fastq_dir  Directory containing paired FASTQ files and
                --pattern    A grep pattern that includes the fastq files
              """
        }

      }

      if (params.prep == "paired-end") {

        if (params.sample_file) {
            Channel
                  .fromPath(params.sample_file, checkIfExists: true)
                  .splitCsv(header: true, strip: true)
                  .map { row ->
                        def sample = row.SampleName
                        def r1     = file(row.Read1, checkIfExists: true)
                        def r2     = file(row.Read2, checkIfExists: true)

                        if (!sample) error "Missing 'SampleName' field in row: ${row}"
                        if (!r1)     error "Read1 not found for sample: ${sample}"
                        if (!r2)     error "Read2 not found for sample: ${sample}"
                        return [sample, tuple(r1, r2)]

                  }
                  .set { raw_reads_ch }

        } else if (!params.sample_file && params.fastq_dir && params.pattern) {
          Channel
                .fromFilePairs("${params.fastq_dir}/${params.pattern}", checkIfExists: true)
                .ifEmpty { exit 1, "You provided an empty directory" }
                .map { file -> tuple(file[0], file[1])}
                .set { raw_reads_ch }
        } else {
            error """
              No input provided. Please specify either:
                --input      Path to a CSV samplesheet OR
                --fastq_dir  Directory containing paired FASTQ files and
                --pattern    A grep pattern that includes the fastq files
              """
        }
      }

      // STAR index construction
      if (!params.genome_index) {

            if (params.genome_fasta_file && params.genome_annot_file) {

                star_index_params_ch = Channel.of([params.genome_fasta_file, params.genome_annot_file])
                STAR_INDEX ( star_index_params_ch, params.tx2gene_extraction_script )
                star_index_ch = STAR_INDEX.out.star_index
                genome_id_file = STAR_INDEX.out.tx2gene_file

            }

      } else if (params.genome_index){

            star_index_ch = params.genome_index

            if (!params.genome_id_file) {
                exit 1, "Please provide a table with corresponding transcript and gene IDs."
            } else {
                genome_id_file = params.genome_id_file
            }

      } else {

            exit 1, """Please provide the necessary files for STAR to construct an index.
            They are:
            --genome_fasta_file
            --genome_annot_file

            OR provide a STAR index --genome_index"""

      }

      // Salmon index construction
      if (!params.salmon_index) {

            if (params.transcriptome_file && params.genome_fasta_file && params.genome_annot_file) {

                SALMON_INDEX ( params.transcriptome_file, params.genome_fasta_file, params.genome_annot_file, params.decoy_script )
                salmon_index_ch = SALMON_INDEX.out.sal_index
            }

      } else if (params.salmon_index) {

            salmon_index_ch = params.salmon_index

      } else {

            exit 1, """Please provide the necessary files for salmon to construct an index.
            They are:
            --transcriptome_file
            --genome_fasta_file
            --genome_annot_file

            OR provide a Salmon index --salmon_index"""

      }

      // Initial QC with fastqc and multiqc
      FASTQC( raw_reads_ch, "initialQC" )
      initial_fastqc_ch = FASTQC.out.fastqc_out.collect()
      MULTIQC( initial_fastqc_ch, params.multiqc_config, "initialQC", "fastqc" )

      // trimming with cutadapt
      CUTADAPT( raw_reads_ch, params.adapters )
      cleanedQC( CUTADAPT.out.trimmed_fastqs, "cleanedQC" )

      // Gathering log for second multiqc
      post_trimming_logs_ch = cleanedQC.out.fastqc_out
        .mix(CUTADAPT.out.trim_logs)
        .collect()
      cleanedMULTIQC( post_trimming_logs_ch, params.multiqc_config, "cleanedQC", "fastqc")

      // Quantification
      // Salmon
      SALMON( CUTADAPT.out.trimmed_fastqs, salmon_index_ch, params.stranded_script )

      salmon_logs_ch = cleanedQC.out.fastqc_out
        .mix(SALMON.out.sal_logs, CUTADAPT.out.trim_logs)
        .collect()

      salmonMULTIQC( salmon_logs_ch, params.multiqc_config, "salmonResults", "salmon" )

      // Mapping against genomes
      if (params.splice_junctions) {

          STAR_FIRST_PASS( CUTADAPT.out.trimmed_fastqs, params.genome_index )
          MERGE_SJ( STAR_FIRST_PASS.out.SJ_logs.map {sid, sj -> sj}.collect(), params.min_unique_reads )
          STAR_SECOND_PASS( CUTADAPT.out.trimmed_fastqs, params.genome_index, MERGE_SJ.out.merged_sj )
          star_map_ch = STAR_SECOND_PASS.out.bams
          star_logs_ch = STAR_SECOND_PASS.out.star_logs

      } else {

          STAR( CUTADAPT.out.trimmed_fastqs, star_index_ch )
          star_map_ch = STAR.out.bams
          star_logs_ch = STAR.out.star_logs

      }


      RIBOSTAR( CUTADAPT.out.trimmed_fastqs, params.rRNA_index, "riboMapping" )
      SPACERSTAR( CUTADAPT.out.trimmed_fastqs, params.rRNA_Spacers_index, "spacerMapping" )

      genome_logs_ch = cleanedQC.out.fastqc_out
        .mix( CUTADAPT.out.trim_logs, star_logs_ch )
        .collect()

      ribo_logs_ch = cleanedQC.out.fastqc_out
        .mix( CUTADAPT.out.trim_logs, RIBOSTAR.out.star_logs )
        .collect()

      spacer_logs_ch = cleanedQC.out.fastqc_out
        .mix( CUTADAPT.out.trim_logs, SPACERSTAR.out.star_logs )
        .collect()

      // Mapping QC
      starMULTIQC( genome_logs_ch, params.multiqc_config, "genomeMultiQC", "star" )
      riboMULTIQC( ribo_logs_ch, params.multiqc_config, "riboMULTIQC", "star" )
      spacerMULTIQC( spacer_logs_ch, params.multiqc_config, "spacerMULTIQC", "star" )

      // Star quantification - running on exons and genes, stranded and unstranded
      READCOUNTS( star_map_ch, params.genome_annot_file, SALMON.out.strand_value )


      READCOUNTS.out.FC_summary_files
          .collect().toSortedList().flatten()
          .branch { file ->
            exonStranded: file.name.contains('exon.stranded')
            exonUnstranded: file.name.contains('exon.unstranded')
            geneStranded: file.name.contains('gene.stranded')
            geneUnstranded: file.name.contains('gene.unstranded')
          }
          .set { FC_files }


      exonStranded_ch = FC_files.exonStranded
          .collect()
          .map { files -> tuple("exonStranded", files) }

      exonUnstranded_ch = FC_files.exonUnstranded
          .collect()
          .map { files -> tuple("exonUnstranded", files) }

      geneStranded_ch = FC_files.geneStranded
          .collect()
          .map { files -> tuple("geneStranded", files) }

      geneUnstranded_ch = FC_files.geneUnstranded
          .collect()
          .map { files -> tuple("geneUnstranded", files) }

      combined_channel = exonStranded_ch.mix(exonUnstranded_ch, geneStranded_ch, geneUnstranded_ch)

      // Quantification QC
      FEATURECOUNTSMULTIQC( combined_channel, params.multiqc_config )

      all_logs_ch = cleanedQC.out.fastqc_out
        .mix( CUTADAPT.out.trim_logs, star_logs_ch, SALMON.out.sal_logs)
        .collect()


      // Gather all logs together
      finalMULTIQC_STRANDED(exonStranded_ch.mix(geneStranded_ch), all_logs_ch, params.multiqc_config )
      finalMULTIQC_UNSTRANDED(exonUnstranded_ch.mix(geneUnstranded_ch), all_logs_ch, params.multiqc_config )

      // Combine all salmon quant files
      salmon_quant_ch = SALMON.out.quant_file.collect()
      GATHERGENEIDS( genome_id_file, params.tx2gene_salmon_file, salmon_quant_ch )

      // Print out all multiqc txt files
      total_multiqc_ch = MULTIQC.out.initial_log
          .mix( cleanedMULTIQC.out.cleaned_log, cleanedMULTIQC.out.cutadapt_log, salmonMULTIQC.out.salmon_log, starMULTIQC.out.genome_log, riboMULTIQC.out.ribo_log, spacerMULTIQC.out.spacer_log )
          .mix( FEATURECOUNTSMULTIQC.out.exonStranded_out, FEATURECOUNTSMULTIQC.out.geneStranded_out, FEATURECOUNTSMULTIQC.out.exonUnstranded_out, FEATURECOUNTSMULTIQC.out.geneUnstranded_out )
          .collect()

      GATHERLOGS( total_multiqc_ch )

      // Gather all multiqc txt files to merge into 1 summary table
      MERGEMULTIQCS( total_multiqc_ch, params.merger_file )


}
