library(tximport)
library(tximportData)
library(data.table)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript your_script.R <quant_dir> <tx2gene_file>", call. = FALSE)
}

input_dir <- args[1]

tx2gene_file <- args[2]

quant_files <- list.files(input_dir, ".sf")

sample_names <- sapply(quant_files, function(x){
  new_name <- gsub("_quant.sf", "", x)
  return(new_name)
})

names(quant_files) <- sample_names

tx2gene_tab <- fread(tx2gene_file, sep = "\t", data.table = FALSE)

path2quant <- sapply(quant_files, function(x){
  paste0(input_dir, "/", x)
})

txi <- tximport(path2quant, type = "salmon", tx2gene = tx2gene_tab, txIdCol = "Transcript stable ID version", geneIdCol = "Gene stable ID")

gene.exp.table <- as.data.frame(txi$counts)

final.table <- merge(tx2gene_tab, gene.exp.table, by.x = "Gene stable ID", by.y = 0, all.y = TRUE)

final.table <- as.data.frame(final.table[!duplicated(final.table$`Gene stable ID`),])

write.table(final.table, "gene_expressions_values.tsv", row.names = FALSE, sep = "\t", quote = FALSE)
