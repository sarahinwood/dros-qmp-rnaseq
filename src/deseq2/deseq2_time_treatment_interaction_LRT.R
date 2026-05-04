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

## power analysis on dds containing samples from both batches except for 120h where only batch 2 is used (batch 1 sample clustering is inconsistent)

dds <- readRDS(dds_file)
gtf <- readGFF(gtf_file)
gtf_gene <- subset(gtf, type=="gene")
gene_to_name <- gtf_gene[,c(9, 10)]

# subset samples that are less comparable
dds <- dds[,dds$treatment_type=="standard"]

## set factors
dds$treatment <- factor(dds$treatment, levels = c("control", "QMP"))
dds$batch <- factor(dds$batch)
dds$timepoint <- factor(dds$timepoint, levels = c("12", "24", "48", "120"))
dds_coldata <- data.frame(colData(dds))
dds_res <- dds

# run DESeq
design(dds_res) <- ~ batch + timepoint + treatment + timepoint:treatment
# reduced model will remove:
    # the effect of time
    # and the effect of QMP which does not vary over time
dds_res <- DESeq(dds_res, test="LRT", reduced = ~ batch + timepoint + treatment)
saveRDS(dds_res, snakemake@output[["dds"]])

res_group_0.1 <- results(dds_res, alpha = 0.1)
paste("timepoint:treatment interaction with padj<0.1")
summary(res_group_0.1)

res_group_0.05 <- results(dds_res, alpha = 0.05)
paste("timepoint:treatment interaction with padj<0.05")
summary(res_group_0.05)

res_group_0.01 <- results(dds_res, alpha = 0.01)
paste("timepoint:treatment interaction with padj<0.01")
summary(res_group_0.01)

ordered_res_group <- res_group_0.05[order(res_group_0.05$padj),]
##Make data table
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
fwrite(ordered_res_group_table, snakemake@output[["res_interaction"]])
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
sig_annots_interaction <- merge(ordered_sig_res_group_table, gene_to_name, by.x="rn", by.y="gene_id")
fwrite(sig_annots_interaction, snakemake@output[["sig_res_interaction"]])
## these are all genes that have a time-specific effect of QMP treatment

# write log
sessionInfo()
