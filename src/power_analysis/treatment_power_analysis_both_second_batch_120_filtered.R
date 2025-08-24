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

library(RnaSeqSampleSize)
library(DESeq2)
library(tidyverse)
library(data.table)

###########
# GLOBALS #
###########

dds_file <- snakemake@input[["dds_file"]]

########
# MAIN #
########

dds <- readRDS(dds_file)
counts <- data.frame(counts(dds, normalized=F))
paste("dds read in, counts df generated")

# make table to fill with res - testing different LFC thresholds
power_res <- data.frame(matrix(ncol = 4, nrow = 1))
col_names <- c("LFC1", "LFC2", "LFC5", "LFC10")
colnames(power_res) <- col_names

## treatment
# reorder so table ordered by factor being tested - etoh vs qmp (not including qmp then etoh samples)
treatment_counts <- counts[,c(
                            ## ETOH samples
                            1:4,9:10,13:14,19:20,
                            ## QMP samples
                            5:8,11:12,17:18,21:22)]
paste("treatment count table generated")

treatment_res <- power_res
# indicate the two sample groups (should be ordered in table)
treatment_disp <- est_count_dispersion(treatment_counts, group=c(rep(0,10), rep(1,10)))
paste("est_count_dispersion finished")

## estimate power with different LFC thresholds

# n = sample number - from online documentation I think this is sample number in each group not total sample number!
# f = alpha threshold
# rho - LFC threshold
# repNumber = number of genes used in estimation

treatment_res$LFC1 <- est_power_distribution(n=10, f=0.05, rho=1, distributionObject=treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC1 complete")
treatment_res$LFC2 <- est_power_distribution(n=10, f=0.05, rho=2, distributionObject=treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC2 complete")
treatment_res$LFC5 <- est_power_distribution(n=10, f=0.05, rho=5, distributionObject=treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC5 complete")
treatment_res$LFC10 <- est_power_distribution(n=10, f=0.05, rho=10, distributionObject=treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC10 complete")
rownames(treatment_res) <- "treatment"
treatment_res <- data.table(treatment_res, keep.rownames=T)

# adjust values from proportion to percentage
treatment_res$LFC1 <- 100*treatment_res$LFC1
treatment_res$LFC2 <- 100*treatment_res$LFC2
treatment_res$LFC5 <- 100*treatment_res$LFC5
treatment_res$LFC10 <- 100*treatment_res$LFC10

fwrite(treatment_res, snakemake@output[["power_res"]])

# write log
sessionInfo()