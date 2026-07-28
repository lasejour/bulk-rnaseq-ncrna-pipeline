#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Please provide the location of the folder where the sample was run", call.=FALSE)
} else if (length(args)==1) {
  #print(args[1])
  setwd(args[1])
}

ReQ_packages = c("dplyr", "data.table")

for (pack in 1:length(ReQ_packages)) {
  if(ReQ_packages[pack] %in% rownames(installed.packages()) == FALSE) {
    install.packages(ReQ_packages[pack], repos='http://cran.us.r-project.org')
  }
}

library(data.table)
library(dplyr)

getwd()

files<-c('initial_multiqc_fastqc.txt', 'cleaned_multiqc_fastqc.txt', 'multiqc_cutadapt.txt', 'genome_multiqc_star.txt', 'ribo_multiqc_star.txt', 'spacer_multiqc_star.txt','multiqc_salmon.txt', 'exonStranded_multiqc_featureCounts.txt', 'exonUnstranded_multiqc_featureCounts.txt', 'geneStranded_multiqc_featureCounts.txt','geneUnstranded_multiqc_featureCounts.txt') 

for (f in files) {
  if (file.exists(f) !=TRUE) {
    print(paste(f,'does not exist. Merging the rest.'))
  } else {
    assign(gsub('.txt','',f), fread(f))
  }
}

# Parse out file extensions from STAR logs
genome_multiqc_star$Sample <- sub('.Log*', '', genome_multiqc_star$Sample, perl = TRUE) 
ribo_multiqc_star$Sample <- sub('-ribo-.*', '', ribo_multiqc_star$Sample, perl = TRUE)
spacer_multiqc_star$Sample <- sub('-spacer-.*', '', spacer_multiqc_star$Sample, perl = TRUE)

# 
colnames(multiqc_cutadapt)[1]<-"cutadapt_Sample"
merged_tables <- cbind(initial_multiqc_fastqc, multiqc_cutadapt[match(initial_multiqc_fastqc$Sample, multiqc_cutadapt$cutadapt_Sample)])
merged_tables$cutadapt_Sample<-NULL

cleaned_multiqc_fastqc$`File type` <- NULL
cleaned_multiqc_fastqc$Encoding <- NULL
colnames(cleaned_multiqc_fastqc) <- paste0('cleaned_',colnames(cleaned_multiqc_fastqc))
merged_tables <- cbind(merged_tables, cleaned_multiqc_fastqc[match(merged_tables$Sample, cleaned_multiqc_fastqc$cleaned_Sample)])
merged_tables$cleaned_Sample<-NULL

colnames(genome_multiqc_star) <- paste0('genome_',colnames(genome_multiqc_star))
merged_tables <- cbind(merged_tables, genome_multiqc_star[match(merged_tables$Sample, genome_multiqc_star$genome_Sample)])
merged_tables$genome_Sample <- NULL

colnames(ribo_multiqc_star) <- paste0('ribo_',colnames(ribo_multiqc_star))
merged_tables <- cbind(merged_tables, ribo_multiqc_star[match(merged_tables$Sample, ribo_multiqc_star$ribo_Sample)])
merged_tables$ribo_Sample <- NULL

colnames(spacer_multiqc_star) <- paste0('spacer_',colnames(spacer_multiqc_star))
merged_tables <- cbind(merged_tables, spacer_multiqc_star[match(merged_tables$Sample, spacer_multiqc_star$spacer_Sample)])
merged_tables$spacer_Sample <- NULL

colnames(multiqc_salmon) <- paste0('salmon_',colnames(multiqc_salmon))
colnames(multiqc_salmon)[2]<- 'salmon_version'
merged_tables <- cbind(merged_tables, multiqc_salmon[match(merged_tables$Sample, multiqc_salmon$salmon_Sample)])
merged_tables$salmon_Sample <- NULL
merged_tables$salmon_start_time <- NULL
merged_tables$salmon_end_time <- NULL
merged_tables$salmon_index_name_hash <- NULL
merged_tables$salmon_index_seq_hash <- NULL
merged_tables$salmon_index_seq_hash512 <- NULL
merged_tables$salmon_index_name_hash512 <- NULL

colnames(exonStranded_multiqc_featureCounts) <- paste0('exon_str_',colnames(exonStranded_multiqc_featureCounts))
exonStranded_multiqc_featureCounts$exon_str_Sample <- sub('_R1_001', '', exonStranded_multiqc_featureCounts$exon_str_Sample, perl = TRUE)
merged_tables <- cbind(merged_tables, exonStranded_multiqc_featureCounts[match(merged_tables$Sample, exonStranded_multiqc_featureCounts$exon_str_Sample)])
merged_tables$exon_str_Sample <- NULL

colnames(exonUnstranded_multiqc_featureCounts) <- paste0('exon_unstr_',colnames(exonUnstranded_multiqc_featureCounts))
exonUnstranded_multiqc_featureCounts$exon_unstr_Sample <- sub('_R1_001', '', exonUnstranded_multiqc_featureCounts$exon_unstr_Sample, perl = TRUE)
merged_tables <- cbind(merged_tables, exonUnstranded_multiqc_featureCounts[match(merged_tables$Sample, exonUnstranded_multiqc_featureCounts$exon_unstr_Sample)])
merged_tables$exon_unstr_Sample <- NULL

colnames(geneStranded_multiqc_featureCounts) <- paste0('gene_str_',colnames(geneStranded_multiqc_featureCounts))
geneStranded_multiqc_featureCounts$gene_str_Sample <- sub('_R1_001', '', geneStranded_multiqc_featureCounts$gene_str_Sample, perl = TRUE)
merged_tables <- cbind(merged_tables, geneStranded_multiqc_featureCounts[match(merged_tables$Sample, geneStranded_multiqc_featureCounts$gene_str_Sample)])
merged_tables$gene_str_Sample <- NULL

colnames(geneUnstranded_multiqc_featureCounts) <- paste0('gene_unstr_',colnames(geneUnstranded_multiqc_featureCounts))
geneUnstranded_multiqc_featureCounts$gene_unstr_Sample <- sub('_R1_001', '', geneUnstranded_multiqc_featureCounts$gene_unstr_Sample, perl = TRUE)
merged_tables <- cbind(merged_tables, geneUnstranded_multiqc_featureCounts[match(merged_tables$Sample, geneUnstranded_multiqc_featureCounts$gene_unstr_Sample)])
merged_tables$gene_unstr_Sample <- NULL

# # Get Metrics
merged_tables %<>% as_tibble() %>% mutate(.,
              Gene_Str_Unstr_Ratio = gene_str_Assigned/gene_unstr_Assigned,
              Exon_Str_Unstr_Ratio = exon_str_Assigned/exon_unstr_Assigned,
              
              Percent_Ribosomal_vs_initial = 100 *(ribo_uniquely_mapped+ribo_multimapped)/`Total Sequences`,
              Percent_Ribosomal_vs_trimmed = 100 *(ribo_uniquely_mapped+ribo_multimapped)/`cleaned_Total Sequences`,
          
              # 
              Percent_Genes_Stranded_uniquely_mapped = 100 * (gene_str_Assigned/genome_uniquely_mapped),
              Percent_Genes_Unstranded_uniquely_mapped = 100 * (gene_unstr_Assigned/genome_uniquely_mapped),
              Percent_Exon_Stranded_uniquely_mapped = 100 * (exon_str_Assigned/genome_uniquely_mapped),
              Percent_Exon_Unstranded_uniquely_mapped = 100 * (exon_unstr_Assigned/genome_uniquely_mapped),
              # 
              Intergenic_Reads_Stranded_uniquely_mapped = genome_uniquely_mapped - gene_str_Assigned,
              Introns_Reads_Stranded_uniquely_mapped = gene_str_Assigned - exon_str_Assigned,
              # 
              Intergenic_Reads_Unstranded_uniquely_mapped = genome_uniquely_mapped - gene_unstr_Assigned,
              Introns_Reads_Unstranded_uniquely_mapped = gene_unstr_Assigned - exon_unstr_Assigned,
              # 
              Percent_Intergenic_Stranded_uniquely_mapped = 100 * (Intergenic_Reads_Stranded_uniquely_mapped/genome_uniquely_mapped),
              Percent_Exonic_Stranded_uniquely_mapped = 100 * (exon_str_Assigned/genome_uniquely_mapped),
              Percent_Intronic_Stranded_uniquely_mapped = 100 * (Introns_Reads_Stranded_uniquely_mapped/genome_uniquely_mapped),
              # 
              Percent_Intergenic_Unstranded_uniquely_mapped = 100 * (Intergenic_Reads_Unstranded_uniquely_mapped/genome_uniquely_mapped),
              Percent_Exonic_Unstranded_uniquely_mapped = 100 * (exon_unstr_Assigned/genome_uniquely_mapped),
              Percent_Intronic_Unstranded_uniquely_mapped = 100 * (Introns_Reads_Unstranded_uniquely_mapped/genome_uniquely_mapped),
              # 
              # 
              # # #######################################################
              # 
              Percent_Genes_Stranded_total_sequences = 100 * (gene_str_Assigned/`Total Sequences`),
              Percent_Genes_Unstranded_total_sequences = 100 * (gene_unstr_Assigned/`Total Sequences`),

              Percent_Exon_Stranded_total_sequences = 100 * (exon_str_Assigned/`Total Sequences`),
              Percent_Exon_Unstranded_total_sequences = 100 * (exon_unstr_Assigned/`Total Sequences`),

              Intergenic_Reads_Stranded_total_sequences = `Total Sequences` - gene_str_Assigned,
              Introns_Reads_Stranded_total_sequences = gene_str_Assigned - exon_str_Assigned,

              Intergenic_Reads_Unstranded_total_sequences = `Total Sequences` - gene_unstr_Assigned,
              Introns_Reads_Unstranded_total_sequences = gene_unstr_Assigned - exon_unstr_Assigned,

              Percent_Intergenic_Stranded_total_sequences = 100 * (Intergenic_Reads_Stranded_total_sequences/`Total Sequences`),
              Percent_Exonic_Stranded_total_sequences = 100 * (exon_str_Assigned/`Total Sequences`),
              Percent_Intronic_Stranded_total_sequences = 100 * (Introns_Reads_Stranded_total_sequences/`Total Sequences`),

              Percent_Intergenic_Unstranded_total_sequences = 100 * (Intergenic_Reads_Unstranded_total_sequences/`Total Sequences`),
              Percent_Exonic_Unstranded_total_sequences = 100 * (exon_unstr_Assigned/`Total Sequences`),
              Percent_Intronic_Unstranded_total_sequences = 100 * (Introns_Reads_Unstranded_total_sequences/`Total Sequences`),

              #########################################################

              Percent_Genes_Stranded_cleaned_sequences = 100 * (gene_str_Assigned/`cleaned_Total Sequences`),
              Percent_Genes_Unstranded_cleaned_sequences = 100 * (gene_unstr_Assigned/`cleaned_Total Sequences`),

              Percent_Exon_Stranded_cleaned_sequences = 100 * (exon_str_Assigned/`cleaned_Total Sequences`),
              Percent_Exon_Unstranded_cleaned_sequences = 100 * (exon_unstr_Assigned/`cleaned_Total Sequences`),

              Intergenic_Reads_Stranded_cleaned_sequences = `cleaned_Total Sequences` - gene_str_Assigned,
              Introns_Reads_Stranded_cleaned_sequences = gene_str_Assigned - exon_str_Assigned,

              Intergenic_Reads_Unstranded_cleaned_sequences = `cleaned_Total Sequences` - gene_unstr_Assigned,
              Introns_Reads_Unstranded_cleaned_sequences = gene_unstr_Assigned - exon_unstr_Assigned,

              Percent_Intergenic_Stranded_cleaned_sequences = 100 * (Intergenic_Reads_Stranded_cleaned_sequences/`cleaned_Total Sequences`),
              Percent_Exonic_Stranded_cleaned_sequences = 100 * (exon_str_Assigned/`cleaned_Total Sequences`),
              Percent_Intronic_Stranded_cleaned_sequences = 100 * (Introns_Reads_Stranded_cleaned_sequences/`cleaned_Total Sequences`),

              Percent_Intergenic_Unstranded_cleaned_sequences = 100 * (Intergenic_Reads_Unstranded_cleaned_sequences/`cleaned_Total Sequences`),
              Percent_Exonic_Unstranded_cleaned_sequences = 100 * (exon_unstr_Assigned/`cleaned_Total Sequences`),
              Percent_Intronic_Unstranded_cleaned_sequences = 100 * (Introns_Reads_Unstranded_cleaned_sequences/`cleaned_Total Sequences`)


              )

for (f in files) {
  assign(gsub('.txt','',f), NULL)
}

write.table(merged_tables,file='merged_results.tsv', sep = "\t", row.names = FALSE)
