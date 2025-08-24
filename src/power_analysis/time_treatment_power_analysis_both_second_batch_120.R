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

## testing each time-point treatment pair as pairwise (can't test power with an LRT test)

###########
## 12 hr ##
###########

paste("12hr power test")

hr12_counts <- counts[,c( ## 12hr control
                          1:4,
                          ## 12hr treatment
                          5:8)]

hr12_treatment_res <- power_res
# indicate the two sample groups (should be ordered in table)
hr12_treatment_disp <- est_count_dispersion(hr12_counts, group=c(rep(0,4), rep(1,4)))

## estimate power with different LFC thresholds

# n = sample number - from online documentation I think this is sample number in each group not total sample number!
# f = alpha threshold
# rho - LFC threshold
# repNumber = number of genes used in estimation

hr12_treatment_res$LFC1 <- est_power_distribution(n=4, f=0.05, rho=1, distributionObject=hr12_treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC1 complete")
hr12_treatment_res$LFC2 <- est_power_distribution(n=4, f=0.05, rho=2, distributionObject=hr12_treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC2 complete")
hr12_treatment_res$LFC5 <- est_power_distribution(n=4, f=0.05, rho=5, distributionObject=hr12_treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC5 complete")
hr12_treatment_res$LFC10 <- est_power_distribution(n=4, f=0.05, rho=10, distributionObject=hr12_treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC10 complete")
rownames(hr12_treatment_res) <- "12hr treatment pairwise"
hr12_treatment_res <- data.table(hr12_treatment_res, keep.rownames=T)

# adjust values from proportion to percentage
hr12_treatment_res$LFC1 <- 100*hr12_treatment_res$LFC1
hr12_treatment_res$LFC2 <- 100*hr12_treatment_res$LFC2
hr12_treatment_res$LFC5 <- 100*hr12_treatment_res$LFC5
hr12_treatment_res$LFC10 <- 100*hr12_treatment_res$LFC10

###########
## 24 hr ##
###########

paste("24hr power test")

hr24_counts <- counts[,c( ## 24hr control
  9,10,
  ## 24hr treatment
  11,12)]

hr24_treatment_res <- power_res
# indicate the two sample groups (should be ordered in table)
hr24_treatment_disp <- est_count_dispersion(hr24_counts, group=c(rep(0,2), rep(1,2)))

hr24_treatment_res$LFC1 <- est_power_distribution(n=2, f=0.05, rho=1, distributionObject=hr24_treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC1 complete")
hr24_treatment_res$LFC2 <- est_power_distribution(n=2, f=0.05, rho=2, distributionObject=hr24_treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC2 complete")
hr24_treatment_res$LFC5 <- est_power_distribution(n=2, f=0.05, rho=5, distributionObject=hr24_treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC5 complete")
hr24_treatment_res$LFC10 <- est_power_distribution(n=2, f=0.05, rho=10, distributionObject=hr24_treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC10 complete")
rownames(hr24_treatment_res) <- "24hr treatment pairwise"
hr24_treatment_res <- data.table(hr24_treatment_res, keep.rownames=T)

# adjust values from proportion to percentage
hr24_treatment_res$LFC1 <- 100*hr24_treatment_res$LFC1
hr24_treatment_res$LFC2 <- 100*hr24_treatment_res$LFC2
hr24_treatment_res$LFC5 <- 100*hr24_treatment_res$LFC5
hr24_treatment_res$LFC10 <- 100*hr24_treatment_res$LFC10

###########
## 48 hr ##
###########

paste("48hr power test")

hr48_counts <- counts[,c( ## 48hr control
  19,20,
  ## 48hr treatment
  21,22)]

hr48_treatment_res <- power_res
# indicate the two sample groups (should be ordered in table)
hr48_treatment_disp <- est_count_dispersion(hr48_counts, group=c(rep(0,2), rep(1,2)))

hr48_treatment_res$LFC1 <- est_power_distribution(n=2, f=0.05, rho=1, distributionObject=hr48_treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC1 complete")
hr48_treatment_res$LFC2 <- est_power_distribution(n=2, f=0.05, rho=2, distributionObject=hr48_treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC2 complete")
hr48_treatment_res$LFC5 <- est_power_distribution(n=2, f=0.05, rho=5, distributionObject=hr48_treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC5 complete")
hr48_treatment_res$LFC10 <- est_power_distribution(n=2, f=0.05, rho=10, distributionObject=hr48_treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC10 complete")
rownames(hr48_treatment_res) <- "48hr treatment pairwise"
hr48_treatment_res <- data.table(hr48_treatment_res, keep.rownames=T)

# adjust values from proportion to percentage
hr48_treatment_res$LFC1 <- 100*hr48_treatment_res$LFC1
hr48_treatment_res$LFC2 <- 100*hr48_treatment_res$LFC2
hr48_treatment_res$LFC5 <- 100*hr48_treatment_res$LFC5
hr48_treatment_res$LFC10 <- 100*hr48_treatment_res$LFC10

############
## 120 hr ##
############

paste("120hr power test")

hr120_counts <- counts[,c( ## 120hr control
  13,14,
  ## 120hr treatment
  17,18)]

hr120_treatment_res <- power_res
# indicate the two sample groups (should be ordered in table)
hr120_treatment_disp <- est_count_dispersion(hr120_counts, group=c(rep(0,2), rep(1,2)))

hr120_treatment_res$LFC1 <- est_power_distribution(n=2, f=0.05, rho=1, distributionObject=hr120_treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC1 complete")
hr120_treatment_res$LFC2 <- est_power_distribution(n=2, f=0.05, rho=2, distributionObject=hr120_treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC2 complete")
hr120_treatment_res$LFC5 <- est_power_distribution(n=2, f=0.05, rho=5, distributionObject=hr120_treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC5 complete")
hr120_treatment_res$LFC10 <- est_power_distribution(n=2, f=0.05, rho=10, distributionObject=hr120_treatment_disp, repNumber=1000)
paste("est_power_distribution on LFC10 complete")
rownames(hr120_treatment_res) <- "120hr treatment pairwise"
hr120_treatment_res <- data.table(hr120_treatment_res, keep.rownames=T)

# adjust values from proportion to percentage
hr120_treatment_res$LFC1 <- 100*hr120_treatment_res$LFC1
hr120_treatment_res$LFC2 <- 100*hr120_treatment_res$LFC2
hr120_treatment_res$LFC5 <- 100*hr120_treatment_res$LFC5
hr120_treatment_res$LFC10 <- 100*hr120_treatment_res$LFC10

timepoint_pairwise_res <- full_join(hr12_treatment_res,full_join(hr24_treatment_res, full_join(hr48_treatment_res,hr120_treatment_res)))

fwrite(timepoint_pairwise_res, snakemake@output[["power_res"]])

# write log
sessionInfo()