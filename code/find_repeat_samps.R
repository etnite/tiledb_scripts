## Find repeated samples in TileDB-VCF instance
##
## This script is intended to find repeated samples stored in a TileDB-VCF instance.
## By this we mean instance where the same genotype has been genotyped multiple
## times.
## 
## The current convention for naming samples in the database is:
##
##   <library_prep_num>.<project>.<state>.<year>.<genotype>
##
## For instance:
##
##   157048.NORGRAINS.IL.21.IL20-1234
##
## for a particular Illinois line tested as part of the 2021 Norgrains project.
## The library prep number is a serial number assigned by the ERSGGL to all DNA
## isolations for genotyping.
##
## To use this script, first generate a list of the genotypes in the database:
##
##   $ tiledbvcf list --uri </path/to/db> > <sample_list.txt>
##
## Then give it to this script as a positional argument:
##
##   $ Rscript find_repeat_samps.R <sample_list.txt>
##
## The output will consist of just the repeated samples in the input list, and
## will be output to the same path as the input file, except with "_repeats_only"
## inserted before the extension
################################################################################

library(tools)

## Read in command-line argument (path to input file)
args <- commandArgs(trailingOnly = TRUE)
samps_file <- args[1]

## Construct output name
out_path <- file_path_sans_ext(samps_file)
out_path <- paste0(out_path, "_repeats_only.txt")

## Isolate repeated samples from input
samps <- read.table(samps_file, header = FALSE, sep = ".")
samps <- samps[order(samps$V1), ]
dups <- unique(samps$V5[duplicated(samps$V5)])
samps <- samps[samps$V5 %in% dups, ]

## Write out subset
write.table(samps, file = out_path, sep = ".", quote = FALSE, col.names = FALSE, row.names = FALSE)
