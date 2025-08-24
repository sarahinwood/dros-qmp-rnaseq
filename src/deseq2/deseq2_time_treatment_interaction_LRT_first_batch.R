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
design(dds_res) <- ~ timepoint + treatment + timepoint:treatment
# reduced model will remove:
    # the effect of time
    # and the effect of QMP which does not vary over time
dds_res <- DESeq(dds_res, test="LRT", reduced = ~ timepoint + treatment)
saveRDS(dds_res, snakemake@output[["dds"]])

res_group <- results(dds_res, alpha = 0.05)
summary(res_group)
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
sig_annots_interaction <- merge(ordered_sig_res_group_table, gene_to_name, by.x="rn", by.y="gene_id")
fwrite(sig_annots_interaction, snakemake@output[["res_interaction"]])
## these are all genes that have a time-specific effect of QMP treatment


## can then pull out results for individual timepoints as contained below
resultsNames(dds_res)


### 12 hr ###
qmp_12hrs_res <- results(dds_res, test="Wald",
                         lfcThreshold = 1, alpha = 0.05,
                         name = "treatment_QMP_vs_control")
                          ## name is different here as 12 is reference level
summary(qmp_12hrs_res)
##Order based of padj
ordered_res_group_12hr <- qmp_12hrs_res[order(qmp_12hrs_res$padj),]
##Make data table
ordered_res_group_12hr_table <- data.table(data.frame(ordered_res_group_12hr), keep.rownames = TRUE)
ordered_sig_res_group_12hr_table <- subset(ordered_res_group_12hr_table, padj < 0.05)
sig_annots_12hr <- merge(ordered_sig_res_group_12hr_table, gene_to_name, by.x="rn", by.y="gene_id")
fwrite(sig_annots_12hr, snakemake@output[["res_12hr"]])


### 24 hr ###
qmp_24hrs_res <- results(dds_res, test="Wald",
                         lfcThreshold = 1, alpha = 0.05,
                         name = "timepoint24.treatmentQMP")
summary(qmp_24hrs_res)
##Order based of padj
ordered_res_group_24hr <- qmp_24hrs_res[order(qmp_24hrs_res$padj),]
##Make data table
ordered_res_group_24hr_table <- data.table(data.frame(ordered_res_group_24hr), keep.rownames = TRUE)
ordered_sig_res_group_24hr_table <- subset(ordered_res_group_24hr_table, padj < 0.05)
sig_annots_24hr <- merge(ordered_sig_res_group_24hr_table, gene_to_name, by.x="rn", by.y="gene_id")
fwrite(sig_annots_24hr, snakemake@output[["res_24hr"]])


### 48 hr ###
qmp_48hrs_res <- results(dds_res, test="Wald",
                         lfcThreshold = 1, alpha = 0.05,
                         name = "timepoint48.treatmentQMP")
summary(qmp_48hrs_res)
##Order based of padj
ordered_res_group_48hr <- qmp_48hrs_res[order(qmp_48hrs_res$padj),]
##Make data table
ordered_res_group_48hr_table <- data.table(data.frame(ordered_res_group_48hr), keep.rownames = TRUE)
ordered_sig_res_group_48hr_table <- subset(ordered_res_group_48hr_table, padj < 0.05)
sig_annots_48hr <- merge(ordered_sig_res_group_48hr_table, gene_to_name, by.x="rn", by.y="gene_id")
fwrite(sig_annots_48hr, snakemake@output[["res_48hr"]])


### 120 hr ###
qmp_120hrs_res <- results(dds_res, test="Wald",
                         lfcThreshold = 1, alpha = 0.05,
                         name = "timepoint120.treatmentQMP")
summary(qmp_120hrs_res)
##Order based of padj
ordered_res_group_120hr <- qmp_120hrs_res[order(qmp_120hrs_res$padj),]
##Make data table
ordered_res_group_120hr_table <- data.table(data.frame(ordered_res_group_120hr), keep.rownames = TRUE)
ordered_sig_res_group_120hr_table <- subset(ordered_res_group_120hr_table, padj < 0.05)
sig_annots_120hr <- merge(ordered_sig_res_group_120hr_table, gene_to_name, by.x="rn", by.y="gene_id")
fwrite(sig_annots_120hr, snakemake@output[["res_120hr"]])

# write log
sessionInfo()
