#!/usr/bin/Rscript

## Find outliers and lack-of-consensus repeats from replicated GBS samples
##
## The purpose of this script is to take an input VCF file, which contains
## replicated samples, and find which of these replicates do not cluster with
## their brethren. Problematic replicates can take two forms:
##   1) When a subset of repeated lines cluster together, but there are one or
##      more outliers
##   2) When all the replicates appear different from each other (i.e. lack of
##      any consensus genotype)
##
## This script uses the SNPRelate package to calculate an identity-by-state (IBS)
## relationship matrix. Calculating this matrix is orders of magnitude faster than
## calculating a distance matrix for large sets of lines. The IBS matrix is then
## later converted to a distance matrix.
##
## The script assumes that sample names in the input VCF are are formatted:
##
##   <library_prep_num>:<project>:<state>:<year>:<genotype>
##
## For instance:
##
##   157048:NORGRAINS:IL:21:IL20-1234
##
## for a particular Illinois line tested as part of the 2021 Norgrains project.
## The library prep number is a serial number assigned by the ERSGGL to all DNA
## isolations for genotyping.
##
## The user defines a cutoff between 0 and 1 for considering lines as outliers - 
## the closer to 0, the more stringent. For instance, setting this to 0.1 will
## cause the script to determine outliers/lack of consensus based on lines
## containing greater than 10% difference in alleles (simple matching).
##
## Output consists of a dataframe listing replicates to retain, and those that are
## outliers or ones which have a complete lack of consensus.
##
## NOTE: SNPRelate creates a .gds file from the input .vcf, in the directory where
## the user chooses to write the output .csv file. This file remains after the 
## script finishes and must be manually removed.
################################################################################

rm(list = ls())
library(SNPRelate)
library(stringr)
library(clstutils)
options(stringsAsFactors = FALSE)


#### User-Defined Constants ####

vcf_file <- "/home/gbg_lab_admin/Array_60TB/Wheat_GBS/Northern_nurseries_Apr2022/filt_P+NUWWSN_25miss_05maf_10het_5dp/all_regions_samp_filt.vcf.gz"

## Distance threshold to classify outliers (0 to 1; closer to 0 = more stringent)
dist_thresh <- 0.1

## Number of threads for calculating IBS matrix
## Generally the IBS calculation is quite fast, so may only need 1
nthreads <- 2

## Path to write output .csv file
out_csv <- "/home/gbg_lab_admin/Array_60TB/Wheat_GBS/Northern_nurseries_Apr2022/filt_P+NUWWSN_25miss_05maf_10het_5dp/repeats_classification.csv"

## Path to output optional IBS matrix .csv file
## Set to NULL to disable
out_mat <- "/home/gbg_lab_admin/Array_60TB/Wheat_GBS/Northern_nurseries_Apr2022/filt_P+NUWWSN_25miss_05maf_10het_5dp/repeats_IBS_matrix.csv"


#### Executable ####

out_dir <- dirname(out_csv)
dir.create(out_dir, recursive = TRUE)

## Create GDS file from VCF
snpgdsVCF2GDS(vcf_file,
              file.path(out_dir, "temp.gds"),
              method = "biallelic.only")
genofile <- snpgdsOpen(file.path(out_dir, "temp.gds"))

## Calculate IBS matrix; perform heirarchical clustering; get dist. matrix
ibs <- snpgdsIBS(genofile, num.thread = nthreads)
hclust <- snpgdsHCluster(ibs, need.mat = TRUE, hang = 0.01)
dist_mat <- hclust$dist

## Make dataframe of repeated genotypes
geno_df <- data.frame("full_name" = colnames(dist_mat),
                      "genotype" = sub("^.*:", "", colnames(dist_mat)))
dups <- unique(geno_df$genotype[duplicated(geno_df$genotype)])
geno_df <- geno_df[geno_df$genotype %in% dups, ]

## Initialize outliers list
out_list <- vector("list", length = length(dups))
names(out_list) <- dups

## Loop through genotypes
for (i in dups) {
  
  ## Subset distance matrix
  gen_sel <- geno_df$full_name[geno_df$genotype == i]
  sub_dist <- dist_mat[gen_sel, gen_sel]
  
  ## Find outliers
  outlie <- findOutliers(sub_dist, cutoff = dist_thresh)
  
  ## If all but one reps are classified as outliers, then there is a lack of
  ## consensus. Otherwise identify outliers (if any)
  if (sum(outlie) == length(outlie) - 1) {
    out_list[[i]] <- data.frame("full_name" = names(outlie),
                                "classification" = "no_consens")
  } else {
    out_list[[i]] <- data.frame("full_name" = names(outlie),
                                "classification" = "retain")
    out_list[[i]]$classification[outlie] <- "outlier"
  }
}

## Rbind list, merge with genotypes DF and write out
out_df <- do.call("rbind", out_list)
write.csv(out_df, out_csv, row.names = FALSE)

## Optionally write out IBS matrix
if (!is.null(out_mat)) {
  write.csv(dist_mat, file = out_mat)
}

snpgdsClose(genofile)
