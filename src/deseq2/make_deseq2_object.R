#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

#############
# LIBRARIES #
#############

library(data.table)
library(DESeq2)
library(rtracklayer)
library(tximport)

###########
# GLOBALS #
###########

gtf_file <- snakemake@input[["gtf_file"]]
sample_data_file <- snakemake@input[["sample_data_file"]]

########
# MAIN #
########

gtf <- readGFF(gtf_file)

# make tx 2 gene file
gtf_mrna <- subset(gtf, type=="mRNA")
tx2gene <- gtf_mrna[,c(11,9)]

##Find all salmon quant files
quant_files <- list.files(path="output/02_salmon", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##assign names to quant files from folder name
names(quant_files) <- gsub(".*/quant_(.+)/.*", "\\1", quant_files)

##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene)

##Import table describing samples
sample_data <- fread(sample_data_file, header=TRUE)
# corrected sample sorting - file naming made this hard
sample_table_sorted <- sample_data[match(names(quant_files), sample_data$sample_name), ]
# check ordered correctly - should return TRUE
all(names(quant_files)==sample_table_sorted$sample_name)
all(colnames(txi) == rownames(sample_table_sorted))

## make dds
dds <- DESeqDataSetFromTximport(txi, sample_table_sorted, ~1)
saveRDS(dds, snakemake@output[["dds_file"]])

## filter out low counts - removes 4000-5000 genes
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds_filtered <- dds[keep,]
saveRDS(dds_filtered, snakemake@output[["dds_filtered_file"]])



# subset to only first batch of sequencing
dds_first_batch <- dds[,dds$batch=="original"]
saveRDS(dds_first_batch, snakemake@output[["dds_first_batch_file"]])

# filtered only first batch of sequencing
dds_first_batch_filtered <- dds_filtered[,dds_filtered$batch=="original"]
saveRDS(dds_first_batch_filtered, snakemake@output[["dds_first_batch_filtered_file"]])



# subset to keep only second batch 120 hr samples - 1st batch cluster a little oddly
dds_both_second_batch_120 <- dds[, !(dds$batch == "original" & dds$timepoint == "120")]
saveRDS(dds_both_second_batch_120, snakemake@output[["dds_both_second_batch_120"]])

# and the same for the filtered dds
dds_both_second_batch_120_filtered <- dds_filtered[, !(dds_filtered$batch == "original" & dds_filtered$timepoint == "120")]
saveRDS(dds_both_second_batch_120_filtered, snakemake@output[["dds_both_second_batch_120_filtered"]])



# write log
sessionInfo()
