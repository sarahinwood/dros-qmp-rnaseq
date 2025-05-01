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


###########
# GLOBALS #
###########

dds_file <- snakemake@input[["dds_file"]]

########
# MAIN #
########

dds <- readRDS(dds_file)

# subset samples that are less comparable
dds <- dds[,dds$treatment_type=="standard"]
dds <- dds[,dds$timepoint!="120"]

## set factors
dds$treatment <- factor(dds$treatment)
dds$treatment_type <- factor(dds$treatment_type)
dds$batch <- factor(dds$batch)
dds$timepoint <- factor(dds$timepoint, levels = c("12", "24", "48"))


dds_treatment_coldata <- data.frame(colData(dds_treatment))

# run DESeq
design(dds_treatment) <- ~timepoint+treatment
dds_treatment <- DESeq(dds_treatment, minReplicatesForReplace=Inf)
saveRDS(dds_treatment, snakemake@output[["dds_file"]])

res_group <- results(dds_treatment, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
resultsNames(dds_treatment)

##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
