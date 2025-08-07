library(clusterProfiler)
library(GO.db)
library(DESeq2)
library(data.table)
library(ggplot2)
library(viridis)

## time:treatment interaction enrichment with all batch 2 samples and batch 1(excluding 120hr)

go_annots <- fread("data/dmel-r6.63_gene_association.fb", skip = 5)
setnames(go_annots, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16"),
         new=c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO ID", " DB:Reference",
               "Evidence", "With (or) From", "Aspect", "DB_Object_Name", "DB_Object_Synonym",
               "DB_Object_Type", "taxon", "Date", "Assigned_by", "Annotation Extension"))
go_to_gene <- go_annots[,c(5,2)]

go_to_name <- go_annots[,c(5)]
go_to_name$GO_term <- Term(go_to_name$`GO ID`)

degs <- fread("output/03_deseq/time_treatment_interaction_LRT/both_batches/both_second_batch_120_filtered/both_second_batch_120_filtered_interaction_sig_annots.csv")

dds_filtered <- readRDS("output/03_deseq/dds_files/dds_both_second_batch_120_filtered.rds")
background_genes <- data.table(rownames(dds_filtered))
fwrite(background_genes, "new_files_for_peter/background_genes_retained_after_expression_filtering.csv")

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
