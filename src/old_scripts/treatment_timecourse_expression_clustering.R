library(DESeq2)
library(DEGreport)

## dds file after time:treatment LRT has run
dds_file = "output/03_deseq/time_treatment_interaction_LRT/both_batches/both_second_batch_120/both_second_batch_120_dds.rds"
dds_lrt <- readRDS(dds_file)


res_LRT <- results(dds_lrt)
padj.cutoff <- 0.05

# Subset the LRT results to return genes with padj < 0.05
sig_res_LRT <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

# Transform counts for data visualization
rld <- rlog(dds_lrt, blind=TRUE)
# Extract the rlog matrix from the object
rld_mat <- assay(rld)

# Obtain rlog values for those significant genes
cluster_rlog <- rld_mat[sig_res_LRT$gene, ]

meta <- colData(dds_lrt)
# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
clusters <- degPatterns(cluster_rlog, metadata = meta, time = "timepoint", col=NULL)

# degPatterns needs a expression matrix, the design experiment and the column used to group samples.


ma = assay(rlog(dds_lrt))[row.names(res)[1:100],]
res <- degPatterns(ma, design, time = "group")
