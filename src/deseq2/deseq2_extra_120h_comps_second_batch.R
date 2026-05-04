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

if (!requireNamespace("org.Dm.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Dm.eg.db", type = "source", ask = FALSE)}

library(data.table)
library(DESeq2)
library(rtracklayer)
library(org.Dm.eg.db)
library(pheatmap)
library(tidyverse)
library(viridis)

###########
# GLOBALS #
###########

dds_file <- snakemake@input[["dds_file"]]
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
dds$timepoint <- factor(dds$timepoint, levels = c("12", "24", "48", "120"))
dds_res <- dds

# run DESeq
design(dds_res) <- ~ treatment
dds_res <- DESeq(dds_res, test="LRT", reduced = ~ 1)

paste("~treatment LRT with padj<0.1")
res_group_0.1 <- results(dds_res, alpha = 0.1)
summary(res_group_0.1)

paste("~treatment LRT with padj<0.05")
res_group <- results(dds_res, alpha = 0.05)
summary(res_group)

## can then pull out results for individual timepoints as contained below
resultsNames(dds_res)
saveRDS(dds_res, snakemake@output[["dds_res"]])

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
# vector of flybase gene symbols (e.g. your "rn" column)
qmp_nc_symbols <- sig_annots_qmp_nc$gene_symbol
# map symbols to gene names
sig_annots_qmp_nc$gene_name <- mapIds(
  org.Dm.eg.db,
  keys    = qmp_nc_symbols,
  column  = "GENENAME",
  keytype = "SYMBOL",
  multiVals = "first")
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
# vector of flybase gene symbols (e.g. your "rn" column)
qmpremoved_nc_symbols <- sig_annots_removed_nc$gene_symbol
# map symbols to gene names
sig_annots_removed_nc$gene_name <- mapIds(
  org.Dm.eg.db,
  keys    = qmpremoved_nc_symbols,
  column  = "GENENAME",
  keytype = "SYMBOL",
  multiVals = "first")
fwrite(sig_annots_removed_nc, snakemake@output[["res_removed_nc"]])


### QMP vs removed ###
paste("qmp vs qmp removed wald with padj<0.1")
removed_qmp_res_0.1 <- results(dds_res, test="Wald",
                          lfcThreshold = 1, alpha = 0.1,
                          contrast=c("treatment","removed", "QMP"))
summary(removed_qmp_res_0.1)
paste("qmp vs qmp removed wald with padj<0.05")
removed_qmp_res <- results(dds_res, test="Wald",
                      lfcThreshold = 1, alpha = 0.05,
                      contrast=c("treatment","removed", "QMP"))
## name is different here as 12 is reference level
summary(removed_qmp_res)
##Order based of padj
ordered_res_group_removed_qmp <- removed_qmp_res[order(removed_qmp_res$padj),]
##Make data table
ordered_res_group_removed_qmp_table <- data.table(data.frame(ordered_res_group_removed_qmp), keep.rownames = TRUE)
ordered_sig_res_group_removed_qmp_table <- subset(ordered_res_group_removed_qmp_table, padj < 0.05)
sig_annots_removed_qmp <- merge(ordered_sig_res_group_removed_qmp_table, gene_to_name, by.x="rn", by.y="gene_id")
# vector of flybase gene symbols (e.g. your "rn" column)
qmpremoved_qmp_symbols <- sig_annots_removed_qmp$gene_symbol
# map symbols to gene names
sig_annots_removed_qmp$gene_name <- mapIds(
  org.Dm.eg.db,
  keys    = qmpremoved_qmp_symbols,
  column  = "GENENAME",
  keytype = "SYMBOL",
  multiVals = "first")
fwrite(sig_annots_removed_qmp, snakemake@output[["res_removed_qmp"]])

##############
## heatmaps ##
##############

##vst transform
vst <- varianceStabilizingTransformation(dds_res, blind=FALSE)
vst_assay_dt <- data.table(assay(vst), keep.rownames=TRUE)

sample_to_label <- data.table(data.frame(colData(dds_res)[,c("treatment", "sample_name")]))
sample_to_label <- sample_to_label %>% remove_rownames %>% column_to_rownames(var="sample_name")
# fix capitalization for plotting
sample_to_label$treatment <- recode(sample_to_label$treatment,
                                    "control" = "Control",
                                    "removed" = "QMP treatment removed",
                                    "QMP"     = "QMP treatment")
sample_to_label <- rename(sample_to_label, Treatment = treatment)
plot_colours <- list(Treatment = c(Control="#fca50a", `QMP treatment removed`="#dd513a", `QMP treatment`="#932667"))

##plot ordering
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

## QMP vs NC ##

##subset for DEGs
vst_degs_qmp_nc <- subset(vst_assay_dt, rn %in% sig_annots_qmp_nc$rn)
##turn first row back to row name
vst_degs_qmp_nc <- vst_degs_qmp_nc %>% remove_rownames %>% column_to_rownames(var="rn")
##remove qmp removed samples
vst_degs_qmp_nc <- vst_degs_qmp_nc[,c(1,2,5,6)]
## turn annotations into gene name (or symbol when uncharacterized/empty)
gene_label_map <- setNames(
  ifelse(is.na(sig_annots_qmp_nc$gene_name) | 
           sig_annots_qmp_nc$gene_name == "uncharacterized protein" |
           sig_annots_qmp_nc$gene_name == "",
         sig_annots_qmp_nc$gene_symbol,
         sig_annots_qmp_nc$gene_name),
  sig_annots_qmp_nc$rn)
rownames(vst_degs_qmp_nc) <- make.unique(gene_label_map[rownames(vst_degs_qmp_nc)])
# retain only used treatments
present_treatments <- unique(sample_to_label[colnames(vst_degs_qmp_nc), "Treatment"])
plot_colours_qmp_nc <- list(
  Treatment = plot_colours$Treatment[names(plot_colours$Treatment) %in% present_treatments])

pdf(snakemake@output[["qmp_nc_heatmap"]])
pheatmap(vst_degs_qmp_nc, cluster_rows=TRUE, cluster_cols=FALSE, clustering_callback=callback, show_rownames=TRUE,
                   annotation_col=sample_to_label, annotation_colors=plot_colours_qmp_nc, annotation_names_col=FALSE,
                   show_colnames = FALSE, border_color=NA, color=viridis(50))
dev.off()

## QMP removed vs NC ##

##subset for DEGs
vst_degs_removed_nc <- subset(vst_assay_dt, rn %in% sig_annots_removed_nc$rn)
##turn first row back to row name
vst_degs_removed_nc <- vst_degs_removed_nc %>% remove_rownames %>% column_to_rownames(var="rn")
##remove qmp removed samples
vst_degs_removed_nc <- vst_degs_removed_nc[,c(1,2,3,4)]
## turn annotations into gene name (or symbol when uncharacterized/empty)
gene_label_map <- setNames(
  ifelse(is.na(sig_annots_removed_nc$gene_name) | 
           sig_annots_removed_nc$gene_name == "uncharacterized protein" |
           sig_annots_removed_nc$gene_name == "",
         sig_annots_removed_nc$gene_symbol,
         sig_annots_removed_nc$gene_name),
  sig_annots_removed_nc$rn)
rownames(vst_degs_removed_nc) <- make.unique(gene_label_map[rownames(vst_degs_removed_nc)])
# retain only used treatments
present_treatments <- unique(sample_to_label[colnames(vst_degs_removed_nc), "Treatment"])
plot_colours_removed_nc <- list(
  Treatment = plot_colours$Treatment[names(plot_colours$Treatment) %in% present_treatments])

pdf(snakemake@output[["removed_nc_heatmap"]])
pheatmap(vst_degs_removed_nc, cluster_rows=TRUE, cluster_cols=FALSE, clustering_callback=callback, show_rownames=TRUE,
         annotation_col=sample_to_label, annotation_colors=plot_colours_removed_nc, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
dev.off()

## QMP removed vs QMP ##

##subset for DEGs
vst_degs_removed_qmp <- subset(vst_assay_dt, rn %in% sig_annots_removed_qmp$rn)
##turn first row back to row name
vst_degs_removed_qmp <- vst_degs_removed_qmp %>% remove_rownames %>% column_to_rownames(var="rn")
##remove qmp removed samples
vst_degs_removed_qmp <- vst_degs_removed_qmp[,c(3,4,5,6)]
## turn annotations into gene name (or symbol when uncharacterized/empty)
gene_label_map <- setNames(
  ifelse(is.na(sig_annots_removed_qmp$gene_name) | 
           sig_annots_removed_qmp$gene_name == "uncharacterized protein" |
           sig_annots_removed_qmp$gene_name == "",
         sig_annots_removed_qmp$gene_symbol,
         sig_annots_removed_qmp$gene_name),
  sig_annots_removed_qmp$rn)
rownames(vst_degs_removed_qmp) <- make.unique(gene_label_map[rownames(vst_degs_removed_qmp)])
# retain only used treatments
present_treatments <- unique(sample_to_label[colnames(vst_degs_removed_qmp), "Treatment"])
plot_colours_removed_qmp <- list(
  Treatment = plot_colours$Treatment[names(plot_colours$Treatment) %in% present_treatments])

pdf(snakemake@output[["removed_qmp_heatmap"]])
pheatmap(vst_degs_removed_qmp, cluster_rows=TRUE, cluster_cols=FALSE, clustering_callback=callback, show_rownames=TRUE,
         annotation_col=sample_to_label, annotation_colors=plot_colours_removed_qmp, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
dev.off()

# write log
sessionInfo()
