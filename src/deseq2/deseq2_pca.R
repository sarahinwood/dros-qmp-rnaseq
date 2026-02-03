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
library(ggplot2)
library(viridis)

###########
# GLOBALS #
###########

dds_file <- snakemake@input[["dds_file"]]

########
# MAIN #
########

dds <- readRDS(dds_file)

## set factors
dds$batch <- factor(dds$batch)
dds$timepoint <- factor(dds$timepoint, levels = c("12", "24", "48", "120"))
dds$treatment <- factor(dds$treatment)
dds$treatment_type <- factor(dds$treatment_type)
dds$time_treatment <- factor(paste(dds$timepoint, dds$treatment))

##transformation --> PCA
# uses the design formula to calculate
# the within-group variability (if blind=FALSE)
# or the across-all-samples variability (if blind=TRUE).
vst <- varianceStabilizingTransformation(dds, blind=TRUE)

## time ##
pca_plot_time <- plotPCA(vst, intgroup=c("timepoint", "batch"), returnData=TRUE)
percentVar_time <- round(100 * attr(pca_plot_time, "percentVar"))
# for comparison with previous graph
#pca_plot1$time_treatment <- factor(pca_plot1$time_treatment,
#                                   levels=c("120 control", "120 QMP", "120 removed",
#                                            "12 control", "12 QMP",
#                                            "24 control", "24 QMP",
#                                            "48 control", "48 QMP"))

## similar clustering to previous PCA though slightly different
# x axis largely = time, 24 and 48 less clear clustering but also reduced sample numbers

pdf(snakemake@output[["PCA_time"]])
ggplot(pca_plot_time, aes(x=PC1, y=PC2, shape=batch, colour=timepoint))+
  geom_point(size=3)+
  scale_color_viridis(discrete=TRUE)+
  labs(shape="Batch", colour="Timepoint (hours)")+
  xlab(paste("PC1:", percentVar_time[1], "% variance")) + 
  ylab(paste("PC2:", percentVar_time[2], "% variance")) + 
  coord_fixed(ratio=1)+
  theme_bw()
dev.off()

## treatment ##

pca_plot_treatment <- plotPCA(vst, intgroup=c("treatment", "batch"), returnData=TRUE)
percentVar_treatment <- round(100 * attr(pca_plot_treatment, "percentVar"))

## similar clustering to previous PCA though slightly different
# x axis largely = time, 24 and 48 less clear clustering but also reduced sample numbers
pdf(snakemake@output[["PCA_treatment"]])
ggplot(pca_plot_treatment, aes(x=PC1, y=PC2, shape=batch, colour=treatment))+
  geom_point(size=3)+
  scale_color_viridis(discrete=TRUE)+
  labs(shape="Batch", colour="Treatment")+
  xlab(paste("PC1:", percentVar_treatment[1], "% variance")) + 
  ylab(paste("PC2:", percentVar_treatment[2], "% variance")) + 
  coord_fixed(ratio=1)+
  theme_bw()
dev.off()

## timepoint & treatment ##

pca_plot_time_treatment <- plotPCA(vst, intgroup=c("timepoint", "treatment", "batch"), returnData=TRUE)
pca_plot_time_treatment$timepoint_treatment <- paste(pca_plot_time_treatment$timepoint, pca_plot_time_treatment$treatment)
pca_plot_time_treatment$timepoint_treatment <- factor(pca_plot_time_treatment$timepoint_treatment,
                                                      levels=c("12 control", "12 QMP", "24 control", "24 QMP",
                                                               "48 control", "48 QMP", "120 control", "120 QMP", "120 removed"))

percentVar_time_treatment <- round(100 * attr(pca_plot_time_treatment, "percentVar"))

## similar clustering to previous PCA though slightly different
# x axis largely = time, 24 and 48 less clear clustering but also reduced sample numbers
pdf(snakemake@output[["PCA_timepoint_treatment"]])
ggplot(pca_plot_time_treatment, aes(x=PC1, y=PC2, shape=batch, colour=timepoint_treatment))+
  geom_point(size=3)+
  scale_color_viridis(discrete=TRUE)+
  labs(shape="Batch", colour="Timepoint (hours) & Treatment")+
  xlab(paste("PC1:", percentVar_time_treatment[1], "% variance")) + 
  ylab(paste("PC2:", percentVar_time_treatment[2], "% variance")) + 
  coord_fixed(ratio=1)+
  theme_bw()
dev.off()

# write log
sessionInfo()

