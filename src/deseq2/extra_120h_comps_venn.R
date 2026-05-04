library(data.table)
library(ggVennDiagram)
library(ggplot2)

qmp_nc <- fread("output/03_deseq/deseq2_extra_120h_comps/second_batch/both_second_batch_120_filtered/both_second_batch_120_filtered_sig_annots_qmp_nc.csv")
qmpremoved_nc <- fread("output/03_deseq/deseq2_extra_120h_comps/second_batch/both_second_batch_120_filtered/both_second_batch_120_filtered_sig_annots_removed_nc.csv")
qmpremoved_qmp <- fread("output/03_deseq/deseq2_extra_120h_comps/second_batch/both_second_batch_120_filtered/both_second_batch_120_filtered_sig_annots_removed_qmp.csv")

gene_sets <- list("QMP treatment\nvs control"      = qmp_nc$rn,
                  "QMP treatment removed\nvs control"  = qmpremoved_nc$rn,
                  "QMP treatment removed\nvs QMP treatment" = qmpremoved_qmp$rn)

venn <- ggVennDiagram(gene_sets, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  theme(legend.position = "none")

ggsave("output/03_deseq/deseq2_extra_120h_comps/second_batch/both_second_batch_120_filtered/venn_DEG_overlap.pdf", venn)


# get gene names in table too

library(org.Dm.eg.db)
library(AnnotationDbi)

# vector of flybase gene symbols (e.g. your "rn" column)
qmp_nc_symbols <- qmp_nc$gene_symbol
# map symbols to gene names
qmp_nc$gene_name <- mapIds(
  org.Dm.eg.db,
  keys    = qmp_nc_symbols,
  column  = "GENENAME",
  keytype = "SYMBOL",
  multiVals = "first")
fwrite(qmp_nc, "output/03_deseq/deseq2_extra_120h_comps/second_batch/both_second_batch_120_filtered/both_second_batch_120_filtered_sig_annots_qmp_nc.csv")

# vector of flybase gene symbols (e.g. your "rn" column)
qmpremoved_nc_symbols <- qmpremoved_nc$gene_symbol
# map symbols to gene names
qmpremoved_nc$gene_name <- mapIds(
  org.Dm.eg.db,
  keys    = qmpremoved_nc_symbols,
  column  = "GENENAME",
  keytype = "SYMBOL",
  multiVals = "first")
fwrite(qmpremoved_nc, "output/03_deseq/deseq2_extra_120h_comps/second_batch/both_second_batch_120_filtered/both_second_batch_120_filtered_sig_annots_removed_nc.csv")

# vector of flybase gene symbols (e.g. your "rn" column)
qmpremoved_qmp_symbols <- qmpremoved_qmp$gene_symbol
# map symbols to gene names
qmpremoved_qmp$gene_name <- mapIds(
  org.Dm.eg.db,
  keys    = qmpremoved_qmp_symbols,
  column  = "GENENAME",
  keytype = "SYMBOL",
  multiVals = "first")
fwrite(qmpremoved_qmp, "output/03_deseq/deseq2_extra_120h_comps/second_batch/both_second_batch_120_filtered/both_second_batch_120_filtered_sig_annots_removed_qmp.csv")





# for LRT results
lrt <- fread("output/03_deseq/time_treatment_interaction_LRT/both_batches/both_second_batch_120_filtered/both_second_batch_120_filtered_interaction_sig_annots.csv")
lrt_symbols <- lrt$gene_symbol
# map symbols to gene names
lrt$gene_name <- mapIds(
  org.Dm.eg.db,
  keys    = lrt_symbols,
  column  = "GENENAME",
  keytype = "SYMBOL",
  multiVals = "first")
fwrite(lrt, "output/03_deseq/time_treatment_interaction_LRT/both_batches/both_second_batch_120_filtered/both_second_batch_120_filtered_interaction_sig_annots.csv")
