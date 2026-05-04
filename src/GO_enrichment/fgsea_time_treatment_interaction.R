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
library(tidyverse)
library(DESeq2)
library(fgsea)
library(ggplot2)
library(viridis)
library(forcats)

###########
# GLOBALS #
###########

go_annots_file <- snakemake@input[["go_annots_file"]]
go_to_name_file <- snakemake@input[["go_to_name_file"]]
all_res_file <- snakemake@input[["all_res_file"]]
background_genes_file <- snakemake@input[["background_genes_file"]]

set.seed(10)

########
# MAIN #
########


### GO term file --> pathways for FGSEA ###

go_annots <- fread(go_annots_file, skip = 5)
setnames(go_annots, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16"),
         new=c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO ID", " DB:Reference",
               "Evidence", "With (or) From", "Aspect", "DB_Object_Name", "DB_Object_Synonym",
               "DB_Object_Type", "taxon", "Date", "Assigned_by", "Annotation Extension"))

# subset to only background genes
background_genes <- fread(background_genes_file)
go_annots_background <- subset(go_annots, DB_Object_ID %in% background_genes$V1)
# table of GO ID to name
go_to_name <- fread(go_to_name_file)
# table with GO terms and gene names as rownames
pathways_list <- go_annots_background %>%
  group_by(`GO ID`) %>%       # Group by GO ID
  summarise(genes = list(unique(DB_Object_ID))) %>%
  deframe()  
# could do this with only biological process terms if I wanted


### DESeq2 results --> ranks for FGSEA ###

all_res <- fread(all_res_file)
## named list with ranks & gene names as names
ranks <- all_res$stat  # Replace with the actual column name/number
names(ranks) <- all_res$rn
ranks <- sort(ranks, decreasing = TRUE)

### run FGSEA & make plot ###
fgsea_res <- fgsea(pathways_list, ranks)
sorted_fgsea_res <- fgsea_res[order(fgsea_res$padj)]
sorted_fgsea_res_no_na <- sorted_fgsea_res[!is.na(padj)]
sum(sorted_fgsea_res_no_na$padj<0.05)
sig_fgsea_res <- subset(sorted_fgsea_res, padj < 0.05)
sig_fgsea_res_annot <- merge(sig_fgsea_res, go_to_name, by.x="pathway", by.y="GO ID")
##need number of genes in leading edge for lollipop plot
sig_fgsea_res_annot$leadingEdge_size <- str_count(sig_fgsea_res_annot$leadingEdge, "FBgn")
fwrite(sig_fgsea_res_annot, snakemake@output[["fgsea_sig_res"]])

plot_res <- sig_fgsea_res_annot %>%
  mutate(GO_term = fct_reorder(GO_term, NES))

pdf(snakemake@output[["fgsea_sig_res_plot"]], width = 10, height = 6)
ggplot(plot_res, aes(x = NES, y = GO_term)) +
  geom_point(aes(size = leadingEdge_size, color = padj)) +
  scale_colour_viridis()+
  scale_size_continuous(name = "Number of genes with\nGO term in leading edge") +
  labs(
    x = "Normalised enrichment score",
    y = "GO Term",
    colour = "Adjusted p-value") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))
dev.off()

# write log
sessionInfo()