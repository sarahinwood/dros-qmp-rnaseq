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
dds$group <- factor(paste(dds$treatment, dds$timepoint, sep = "_"))
dds_wald <- dds

# run DESeq
design(dds_wald) <- ~ batch + group
dds_wald <- DESeq(dds_wald)
saveRDS(dds_wald, snakemake@output[[dds_wald]])

## can then pull out results for individual timepoints as contained below
resultsNames(dds_wald)

### 12h ###
res_12h <- results(dds_wald,
                   test         = "Wald",
                   contrast     = c("group", "QMP_12", "control_12"),
                   lfcThreshold = 1,
                   alpha        = 0.05)
summary(res_12h)
res_12h_ordered <- res_12h[order(res_12h$padj), ]
res_12h_table <- data.table(data.frame(res_12h_ordered), keep.rownames = TRUE)
res_12h_sig <- subset(res_12h_table, padj < 0.05)
res_12h_sig_annot <- merge(res_12h_sig, gene_to_name, by.x = "rn", by.y = "gene_id")
fwrite(res_12h_sig_annot, snakemake@output[["res_12hr"]])

### 24h ###
res_24h <- results(dds_wald,
                   test         = "Wald",
                   contrast     = c("group", "QMP_24", "control_24"),
                   lfcThreshold = 1,
                   alpha        = 0.05)
summary(res_24h)
res_24h_ordered <- res_24h[order(res_24h$padj), ]
res_24h_table <- data.table(data.frame(res_24h_ordered), keep.rownames = TRUE)
res_24h_sig <- subset(res_24h_table, padj < 0.05)
res_24h_sig_annot <- merge(res_24h_sig, gene_to_name, by.x = "rn", by.y = "gene_id")
fwrite(res_24h_sig_annot, snakemake@output[["res_24hr"]])

### 48h ###
res_48h <- results(dds_wald,
                   test         = "Wald",
                   contrast     = c("group", "QMP_48", "control_48"),
                   lfcThreshold = 1,
                   alpha        = 0.05)
summary(res_48h)
res_48h_ordered <- res_48h[order(res_48h$padj), ]
res_48h_table <- data.table(data.frame(res_48h_ordered), keep.rownames = TRUE)
res_48h_sig <- subset(res_48h_table, padj < 0.05)
res_48h_sig_annot <- merge(res_48h_sig, gene_to_name, by.x = "rn", by.y = "gene_id")
fwrite(res_48h_sig_annot, snakemake@output[["res_48hr"]])

### 120h ###
res_120h <- results(dds_wald,
                    test         = "Wald",
                    contrast     = c("group", "QMP_120", "control_120"),
                    lfcThreshold = 1,
                    alpha        = 0.05)
summary(res_120h)
res_120h_ordered <- res_120h[order(res_120h$padj), ]
res_120h_table <- data.table(data.frame(res_120h_ordered), keep.rownames = TRUE)
res_120h_sig <- subset(res_120h_table, padj < 0.05)
res_120h_sig_annot <- merge(res_120h_sig, gene_to_name, by.x = "rn", by.y = "gene_id")
fwrite(res_120h_sig_annot, snakemake@output[["res_120hr"]])

# write log
sessionInfo()