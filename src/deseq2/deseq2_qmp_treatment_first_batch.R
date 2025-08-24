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

###########
# GLOBALS #
###########

dds_file <- snakemake@input[["dds_file"]]
gtf_file <- snakemake@input[["gtf_file"]]

########
# MAIN #
########

dds <- readRDS(dds_file)
gtf <- readGFF(gtf_file)
gtf_gene <- subset(gtf, type=="gene")
gene_to_name <- gtf_gene[,c(9, 10)]

# subset samples that are less comparable
dds <- dds[,dds$treatment_type=="standard"]

## set factors
dds$treatment <- factor(dds$treatment, levels = c("control", "QMP"))
dds$timepoint <- factor(dds$timepoint, levels = c("12", "24", "48", "120"))
dds_coldata <- data.frame(colData(dds))
dds_res <- dds

# run DESeq
design(dds_res) <- ~ timepoint + treatment
dds_res <- DESeq(dds_res, test="Wald")
saveRDS(dds_res, snakemake@output[["dds"]])

resultsNames(dds_res)
res_group <- results(dds_res, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
sig_annots_treatment <- merge(ordered_sig_res_group_table, gene_to_name, by.x="rn", by.y="gene_id")
fwrite(sig_annots_treatment, snakemake@output[["sig_annots_treatment"]])
## these are all genes that are DE between QMP treatment and solvent control

# write log
sessionInfo()
