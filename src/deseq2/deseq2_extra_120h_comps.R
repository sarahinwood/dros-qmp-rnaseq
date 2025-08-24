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

dds_file <- snakemake@input[["dds_file"]] # "output/03_deseq/dds_files/dds_both_batches.rds"
gtf_file <- snakemake@input[["gtf_file"]] # "data/misc_dmel-r6.63_files/dmel-all-r6.63.gtf" 

########
# MAIN #
########

dds <- readRDS(dds_file)
gtf <- readGFF(gtf_file)
gtf_gene <- subset(gtf, type=="gene")
gene_to_name <- gtf_gene[,c(9, 10)]

# subset to 120 hr samples
dds <- dds[,dds$timepoint=="120"]

## set factors
dds$treatment <- factor(dds$treatment, levels = c("control", "removed", "QMP"))
dds$batch <- factor(dds$batch)
dds$timepoint <- factor(dds$timepoint, levels = c("12", "24", "48", "120"))
dds_coldata <- data.frame(colData(dds))
dds_res <- dds

# run DESeq
design(dds_res) <- ~ batch + treatment
# reduced model will remove:
# the effect of batch
dds_res <- DESeq(dds_res, test="LRT", reduced = ~ batch)

paste("~batch+treatment LRT with padj<0.1")
res_group_0.1 <- results(dds_res, alpha = 0.1)
summary(res_group_0.1)

paste("~batch+treatment LRT with padj<0.05")
res_group <- results(dds_res, alpha = 0.05)
summary(res_group)

## can then pull out results for individual timepoints as contained below
resultsNames(dds_res)


### QMP vs control ###
paste("qmp vs nc wald with padj<0.1")
qmp_nc_res_0.1 <- results(dds_res, test="Wald",
                      lfcThreshold = 1, alpha = 0.1,
                      name = "treatment_QMP_vs_control")
## name is different here as 12 is reference level
summary(qmp_nc_res_0.1)
paste("qmp vs nc wald with padj<0.05")
qmp_nc_res <- results(dds_res, test="Wald",
                         lfcThreshold = 1, alpha = 0.05,
                         name = "treatment_QMP_vs_control")
## name is different here as 12 is reference level
summary(qmp_nc_res)
##Order based of padj
ordered_res_group_qmp_nc <- qmp_nc_res[order(qmp_nc_res$padj),]
##Make data table
ordered_res_group_qmp_nc_table <- data.table(data.frame(ordered_res_group_qmp_nc), keep.rownames = TRUE)
ordered_sig_res_group_qmp_nc_table <- subset(ordered_res_group_qmp_nc_table, padj < 0.05)
sig_annots_qmp_nc <- merge(ordered_sig_res_group_qmp_nc_table, gene_to_name, by.x="rn", by.y="gene_id")
fwrite(sig_annots_qmp_nc, snakemake@output[["res_qmp_nc"]])

### removed vs negative control ###
paste("qmp removed vs nc wald with padj<0.1")
removed_nc_res_0.1 <- results(dds_res, test="Wald",
                          lfcThreshold = 1, alpha = 0.1,
                          name = "treatment_removed_vs_control")
## name is different here as 12 is reference level
summary(removed_nc_res_0.1)
paste("qmp removed vs nc wald with padj<0.05")
removed_nc_res <- results(dds_res, test="Wald",
                      lfcThreshold = 1, alpha = 0.05,
                      name = "treatment_removed_vs_control")
## name is different here as 12 is reference level
summary(removed_nc_res)
##Order based of padj
ordered_res_group_removed_nc <- removed_nc_res[order(removed_nc_res$padj),]
##Make data table
ordered_res_group_removed_nc_table <- data.table(data.frame(ordered_res_group_removed_nc), keep.rownames = TRUE)
ordered_sig_res_group_removed_nc_table <- subset(ordered_res_group_removed_nc_table, padj < 0.05)
sig_annots_removed_nc <- merge(ordered_sig_res_group_removed_nc_table, gene_to_name, by.x="rn", by.y="gene_id")
fwrite(sig_annots_removed_nc, snakemake@output[["res_removed_nc"]])


### QMP vs removed ###
paste("qmp vs qmp removed wald with padj<0.1")
qmp_nc_res_0.1 <- results(dds_res, test="Wald",
                      lfcThreshold = 1, alpha = 0.1,
                      contrast=c("treatment","removed", "QMP"))
summary(qmp_nc_res_0.1)
paste("qmp vs qmp removed wald with padj<0.05")
qmp_nc_res <- results(dds_res, test="Wald",
                      lfcThreshold = 1, alpha = 0.05,
                      contrast=c("treatment","removed", "QMP"))
## name is different here as 12 is reference level
summary(qmp_nc_res)
##Order based of padj
ordered_res_group_qmp_nc <- qmp_nc_res[order(qmp_nc_res$padj),]
##Make data table
ordered_res_group_qmp_nc_table <- data.table(data.frame(ordered_res_group_qmp_nc), keep.rownames = TRUE)
ordered_sig_res_group_qmp_nc_table <- subset(ordered_res_group_qmp_nc_table, padj < 0.05)
sig_annots_qmp_nc <- merge(ordered_sig_res_group_qmp_nc_table, gene_to_name, by.x="rn", by.y="gene_id")
fwrite(sig_annots_qmp_nc, snakemake@output[["res_removed_qmp"]])

# write log
sessionInfo()
