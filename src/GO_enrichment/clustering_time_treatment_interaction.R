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

BiocManager::install('DEGreport', type='source', ask=FALSE)

library(data.table)
library(tidyverse)
library(ggplot2)
library(viridis)
library(DESeq2)
library(DEGreport)
library(clusterProfiler)

###########
# GLOBALS #
###########

sample_table_file <- snakemake@input[["sample_table_file"]]
deg_list_file <- snakemake@input[["deg_list_file"]]
dds_file <- snakemake@input[["dds_file"]]
go_annots_file <- snakemake@input[["go_annots_file"]]
go_to_name_file <- snakemake@input[["go_to_name_file"]]
background_genes_file <- snakemake@input[["background_genes_file"]]

set.seed(10)

########
# MAIN #
########

#### cluster into expression patterns ####

sample_table <- fread(sample_table_file)
deg_list <- fread(deg_list_file)
dds <- readRDS(dds_file)
dds_cd <- data.frame(colData(dds))
subset_dds <- dds[,dds$treatment_type=="standard"]
subset_dds_cd <- data.frame(colData(subset_dds))

# Transform counts for data visualization
rld <- rlog(subset_dds, blind=TRUE)
# Extract the rlog matrix from the object
rld_mat <- assay(rld)

# subset to rlog values for those significant genes
cluster_rlog <- rld_mat[deg_list$rn,]

# format metadata correctly
sample_info <- sample_table[,c(1,3,4,5,6)]
sample_info_good_samples <- subset(sample_info, !(sample_info$batch=="original"&sample_info$timepoint=="120"))
sample_info_good_samples <- subset(sample_info_good_samples, treatment_type=="standard")
# order meta table to match rlog and add rownames
meta <- data.frame(sample_info_good_samples)
meta_sorted <- meta[match(colnames(cluster_rlog), meta$sample_name), ]
rownames(meta_sorted) <- meta_sorted$sample_name
meta_sorted$sample_name <- NULL
meta_sorted$batch <- NULL
meta_sorted$treatment_type <- NULL
meta_sorted$timepoint <- factor(meta_sorted$timepoint, levels=c("12", "24", "48", "120"))
meta_sorted$treatment <- factor(meta_sorted$treatment, levels=c("control", "QMP"))

# check samples are ordered correctly
all(colnames(cluster_rlog) %in% rownames(meta_sorted))
all(colnames(cluster_rlog) == rownames(meta_sorted))

clusters <- degPatterns(cluster_rlog, metadata = meta_sorted,
                        time = "timepoint", col="treatment",
                        reduce=TRUE, minc=50)
# clustering settings that can be altered to change results if wanted:
# minc = by default is 15 (minimum number of genes in a group, can be higher or lower)
# reduce = TRUE/FALSE (remove outliers from clusters) - Not used if consensusCluster is TRUE
# consensusCluster = TRUE (FALSE = default) - different clustering method - TRUE takes a long time and returns 2 clusters

pdf(snakemake@output[["cluster_plot"]], width=20, height=20)
clusters$plot
dev.off()

gene_to_cluster <- clusters$df
fwrite(gene_to_cluster, snakemake@output[["genes_to_clusters"]])

# plot showing the different patterns at different values for clustering cuttree function
clusters$benchmarking

# plot showing how the numbers of clusters and genes changed at different values for clustering cuttree function
clusters$benchmarking_curve

## 3732/3823 or 97% DEGs clustered when reduce=false and minc default
## 3108/3823 of 81% DEGs clustered when reduce=false and minc=50

## 3281/3823 or 85% DEGs clustered when reduce=true and minc default
## 2631/3823 or 68% DEGs clustered when reduce=true and minc=50


#### GO term enrichment in each cluster ####

## set up necessary data
go_annots <- fread(go_annots_file, skip = 5)
setnames(go_annots, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16"),
         new=c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO ID", " DB:Reference",
               "Evidence", "With (or) From", "Aspect", "DB_Object_Name", "DB_Object_Synonym",
               "DB_Object_Type", "taxon", "Date", "Assigned_by", "Annotation Extension"))
go_to_gene <- go_annots[,c(5,2)]
go_to_name <- fread(go_to_name_file)
background_genes <- fread(background_genes_file)

# Split your gene_to_cluster into a list by cluster
cluster_genes <- split(gene_to_cluster$genes, gene_to_cluster$cluster)

# Run enrichment for each cluster
cluster_enrichments <- lapply(names(cluster_genes), function(cl){
  message("Running enrichment for cluster ", cl)
  enricher(cluster_genes[[cl]],
           pvalueCutoff = 0.05,
           TERM2GENE = go_to_gene,
           TERM2NAME = go_to_name,
           universe = background_genes$V1)})

# Give names to the list so you can access them later
names(cluster_enrichments) <- names(cluster_genes)

# Filter out clusters with no enrichment results
cluster_enrichments <- Filter(function(x) nrow(as.data.frame(x)) > 0, cluster_enrichments)

# Now safely rbind
all_results <- rbindlist(lapply(names(cluster_enrichments), function(cl){
  df <- as.data.frame(cluster_enrichments[[cl]])
  df$cluster <- cl
  df
}))



all_results_plot <- all_results %>%
  dplyr::mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text = x)))) %>%
  dplyr::arrange(GeneRatio)

all_results_plot$Description <- factor(all_results_plot$Description,
                                  levels = unique(all_results_plot$Description[order(all_results_plot$GeneRatio)]))
all_results_plot$cluster <- as.numeric(all_results_plot$cluster)
all_results_plot$cluster_name <- paste("Cluster", all_results_plot$cluster)
all_results_plot$cluster_name <- factor(all_results_plot$cluster_name, levels=unique(all_results_plot$cluster_name[order(all_results_plot$cluster)]))

fwrite(all_results_plot, snakemake@output[["go_res"]])

pdf(snakemake@output[["go_plot"]], width=50, height=25)
ggplot(all_results_plot, aes(x = GeneRatio, 
                        y = Description, 
                        size = Count, 
                        color = p.adjust)) +
  geom_point() +
  scale_color_viridis() +
  facet_wrap(~ cluster_name, scales = "free_y") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8)) +
  labs(x = "Gene Ratio (Ratio of DEGs in cluster with GO term)", y = "GO term",
       color = "Adjusted p-value", size = "Number of DEGs in\ncluster with GO term")
dev.off()

## could try with Mfuzz for more robust clustering but can't use on mac

# write log
sessionInfo()