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
library(clusterProfiler)
library(ggplot2)
library(viridis)

###########
# GLOBALS #
###########

go_annots_file <- snakemake@input[["go_annots_file"]]
go_to_name_file <- snakemake@input[["go_to_name_file"]]
degs_file <- snakemake@input[["degs_file"]]
dds_filtered_file <- snakemake@input[["dds_filtered_file"]]

set.seed(10)

########
# MAIN #
########

## time:treatment interaction enrichment with all batch 2 samples and batch 1(excluding 120hr)

go_annots <- fread(go_annots_file, skip = 5)
setnames(go_annots, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16"),
         new=c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO ID", " DB:Reference",
               "Evidence", "With (or) From", "Aspect", "DB_Object_Name", "DB_Object_Synonym",
               "DB_Object_Type", "taxon", "Date", "Assigned_by", "Annotation Extension"))
go_to_gene <- go_annots[,c(5,2)]

go_to_name <- fread(go_to_name_file)

degs <- fread(degs_file)

dds_filtered <- readRDS(dds_filtered_file)
background_genes <- data.table(rownames(dds_filtered))
fwrite(background_genes, snakemake@output[["background"]])

results <- enricher(degs$rn,
                    pvalueCutoff = 0.05,
                    TERM2GENE=go_to_gene,
                    TERM2NAME = go_to_name,
                    universe=background_genes$V1)

results_dt <- data.table(as.data.frame(results), keep.rownames = FALSE)

plot_terms <- results_dt %>%
  dplyr::mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text = x)))) %>%
  dplyr::arrange(GeneRatio) %>%
  dplyr::mutate(Description = factor(Description, levels = Description))

pdf(snakemake@output[["plot"]], width = 16, height = 8)
ggplot(plot_terms, aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_colour_viridis()+
  scale_size_continuous(name = "Number of DEGs with GO term") +
  labs(
    x = "Gene Ratio (Ratio of DEGs with GO term)",
    y = "GO Term",
    colour = "Adjusted p-value") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))
dev.off()

# write log
sessionInfo()