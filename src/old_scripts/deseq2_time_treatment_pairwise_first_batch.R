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
dds_res <- DESeq(dds_res)
saveRDS(dds_res, snakemake@output[["dds"]])

res_group <- results(dds_res, lfcThreshold = 1, alpha = 0.1)
summary(res_group)
resultsNames(dds_res)

# For timepoint 12 - 12 is the reference level for time
qmp_12hrs_res <- results(dds_res, lfcThreshold = 1, alpha = 0.1, name = "treatment_QMP_vs_control")
summary(qmp_12hrs_res)
##Order based of padj
ordered_res_group_12hr <- qmp_12hrs_res[order(qmp_12hrs_res$padj),]
##Make data table
ordered_res_group_12hr_table <- data.table(data.frame(ordered_res_group_12hr), keep.rownames = TRUE)
ordered_sig_res_group_12hr_table <- subset(ordered_res_group_12hr_table, padj < 0.1)
sig_annots_12hr <- merge(ordered_sig_res_group_12hr_table, gene_to_name, by.x="rn", by.y="gene_id")
fwrite(sig_annots_12hr, snakemake@output[["res_12hr"]])
##### previously found 0 DEGs #####

# For timepoint 24
qmp_24hrs_res <- results(dds_res, lfcThreshold = 1, alpha = 0.1, name = "timepoint24.treatmentQMP")
summary(qmp_24hrs_res)
##Order based of padj
ordered_res_group_24hr <- qmp_24hrs_res[order(qmp_24hrs_res$padj),]
##Make data table
ordered_res_group_24hr_table <- data.table(data.frame(ordered_res_group_24hr), keep.rownames = TRUE)
ordered_sig_res_group_24hr_table <- subset(ordered_res_group_24hr_table, padj < 0.1)
sig_annots_24hr <- merge(ordered_sig_res_group_24hr_table, gene_to_name, by.x="rn", by.y="gene_id")
fwrite(sig_annots_24hr, snakemake@output[["res_24hr"]])
##### previously found 17 DEGs, 2 downreg #####
### of note ###
## immunity, oxidative stress, signalling pathways, transmembrane transport (voltage-gated calcium channel activity and jap junction components)
## CG14275 unknown function but responsive to notch
## CG1124 potentially responsive to JH
## Ca-alpha1D, Ca-Ma2d, shakB, rgn - nutrition stuff

# For timepoint 48
qmp_48hrs_res <- results(dds_res, lfcThreshold = 1, alpha = 0.1, name = "timepoint48.treatmentQMP")
summary(qmp_48hrs_res)
##Order based of padj
ordered_res_group_48hr <- qmp_48hrs_res[order(qmp_48hrs_res$padj),]
##Make data table
ordered_res_group_48hr_table <- data.table(data.frame(ordered_res_group_48hr), keep.rownames = TRUE)
ordered_sig_res_group_48hr_table <- subset(ordered_res_group_48hr_table, padj < 0.1)
sig_annots_48hr <- merge(ordered_sig_res_group_48hr_table, gene_to_name, by.x="rn", by.y="gene_id")
fwrite(sig_annots_48hr, snakemake@output[["res_48hr"]])
##### previously found 21, 12 downreg DEGs #####
### of note ###
## many with roles in reproduction- often as components of the vitelline membrane or chorion
## altered expression of genes associated with immunity, signalling pathways, transmembrane transport and oxidative stress
## CG14275 again
## Vmat
## CG16926
## Yp3

# For timepoint 120
qmp_120hrs_res <- results(dds_res, lfcThreshold = 1, alpha = 0.1, name = "timepoint120.treatmentQMP")
summary(qmp_120hrs_res)
##Order based of padj
ordered_res_group_120hr <- qmp_120hrs_res[order(qmp_120hrs_res$padj),]
##Make data table
ordered_res_group_120hr_table <- data.table(data.frame(ordered_res_group_120hr), keep.rownames = TRUE)
ordered_sig_res_group_120hr_table <- subset(ordered_res_group_120hr_table, padj < 0.1)
sig_annots_120hr <- merge(ordered_sig_res_group_120hr_table, gene_to_name, by.x="rn", by.y="gene_id")
fwrite(sig_annots_120hr, snakemake@output[["res_120hr"]])
##### previously found 18, 2 downreg DEGs #####
### of note ###
## Cp16
## Obp99b
## Smyd4-2 and tsh
## Cyp6d2

# write log
sessionInfo()
