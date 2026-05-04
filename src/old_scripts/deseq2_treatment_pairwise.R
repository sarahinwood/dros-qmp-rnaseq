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
dds$batch <- factor(dds$batch)
dds$timepoint <- factor(dds$timepoint, levels = c("12", "24", "48", "120"))
dds_coldata <- data.frame(colData(dds))
dds_res <- dds
design(dds_res) <- ~ batch + treatment

## analyse each timepoint separately ##

## 12 hr
dds_res_12hr <- dds_res[,dds_res$timepoint=="12"]
dds_res_12hr <- DESeq(dds_res_12hr)
res_group_12hr <- results(dds_res_12hr, lfcThreshold = 1, alpha = 0.1)
summary(res_group_12hr)

## 24 hr
dds_res_24hr <- dds_res[,dds_res$timepoint=="24"]
design(dds_res_24hr) <- ~treatment
dds_res_24hr <- DESeq(dds_res_24hr)
res_group_24hr <- results(dds_res_24hr, lfcThreshold = 1, alpha = 0.1)
summary(res_group_24hr)

## 48 hr
dds_res_48hr <- dds_res[,dds_res$timepoint=="48"]
design(dds_res_48hr) <- ~treatment
dds_res_48hr <- DESeq(dds_res_48hr)
res_group_48hr <- results(dds_res_48hr, lfcThreshold = 1, alpha = 0.1)
summary(res_group_48hr)

## 120 hr
dds_res_120hr <- dds_res[,dds_res$timepoint=="120"]
dds_res_120hr <- DESeq(dds_res_120hr)
res_group_120hr <- results(dds_res_120hr, lfcThreshold = 1, alpha = 0.1)
summary(res_group_120hr)

# write log
sessionInfo()
