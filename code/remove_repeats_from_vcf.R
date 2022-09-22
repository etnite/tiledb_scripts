## Remove repeats from TileDB-exported VCF file
##
## In the current NorGrains workflow, genotyped samples are given names using the
## following convention:
##
##    <library_prep_num>.<project>.<state>.<year>.<genotype>
##
## For instance:
##
##   157048.NORGRAINS.IL.21.IL20-1234
##
## Therefore genotypes may be repeated within a VCF file if, for instance, they
## were sampled in multiple years. This may or may not be desirable depending on
## what we intend to do with the data. I typically first use the script 
## classify_repeat_samples.R to ensure that I'm first removing low-quality repeated
## samples. After this, there may be cases where we still wish to remove
## all but one repeated genotype.
##
## This script is designed to remove all duplicated genotypes by preferentially
## retaining the duplicate with the lowest amount of missing data for each set of
## repeated genotypes. We can optionally supply the output of classify_repeat_samples.R,
## a file named repeats_classification.csv, to ensure we aren't unwittingly retaining
## a poor-quality replicate.
##
## Finally, we can choose to either retain full sample names in the output, or
## else just the genotype name (i.e. the part after the last colon)
##
## NOTE: bcftools must be present in the user's PATH for this script to work
################################################################################

rm(list = ls())


##### User-Defined Constants ##########

## Paths to input VCF/BCF files
in_vcf_file  <- "/Users/ward.1660/Downloads/OSU_FAES_transfer/remove_repeats_test/all_regions.bcf"
out_vcf_file <- "/Users/ward.1660/Downloads/OSU_FAES_transfer/remove_repeats_test/output/test_out.bcf"

## Path to .csv file created by classify_repeat_samples.R
## Set to NULL to disable checking for bad repeated samples
#repeats_file <- "/Users/ward.1660/Downloads/OSU_FAES_transfer/remove_repeats_test/repeats_classification.csv"
repeats_file <- NULL

## Either "full" to specify that the full sample IDs should be used in the output,
## or "short" to only use the genotype name (after the last colon in the full sample ID)
out_names <- "full"


##### Executable ######################

out_dir <- dirname(out_vcf_file)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Determine if output should be VCF or BCF
if (grepl("bcf", out_vcf_file)) {
  out_fmt <- "b"
  tmp_var_file <- file.path(out_dir, "temp.bcf")
} else if (grepl("vcf", out_vcf_file)) {
  out_fmt <- "v"
  tmp_var_file <- file.path(out_dir, "temp.vcf.gz")
} else {
  stop("Please supply a .bcf or .vcf.gz file for output")
}

## Sanity check - output full or short sample names?
if (!(out_names %in% c("full", "short"))) {
  stop("please select either 'full' or 'short' for out_names constant")
}

## If not using a repeats file, just create the per-sample-counts (PSC file)
## using bcftools stats
if (is.null(repeats_file)) {
  system(sprintf("bcftools stats %s -s - |
                 grep PSC |
                 tail -n +2 > %s",
                 in_vcf_file,
                 file.path(out_dir, "PSC.txt")))
  
## Otherwise remove bad repeat samples before creating PSC file
} else {
  reps <- read.csv(repeats_file)
  drops <- reps$full_name[reps$classification %in% c("no_consens", "outlier")]
  write(drops, file = file.path(out_dir, "drop_samps.txt"), sep = "\n")

  system(sprintf("bcftools view %s -S ^%s --force-samples -Ou |
                 bcftools stats -s - |
                 grep PSC |
                 tail -n +2 > %s",
                 in_vcf_file,
                 file.path(out_dir, "drop_samps.txt"),
                 file.path(out_dir, "PSC.txt")))
  file.remove(file.path(out_dir, "drop_samps.txt"))
}

## Parse the PSC file - order it by missing data content, then eliminate all
## repeated genotypes
psc <- read.table(file.path(out_dir, "PSC.txt"), sep = "\t", header = TRUE, 
                  comment.char = "")
psc <- psc[order(psc$X.14.nMissing), ]
psc$geno <- sub("^.*\\.", "", psc$X.3.sample)
psc <- psc[!duplicated(psc$geno), ]
psc <- psc[order(psc$X.3.sample), ]
write(psc$X.3.sample, file = file.path(out_dir, "samps_subset.txt"), sep = "\n")
write.table(psc[, c("X.3.sample", "geno")],
            file = file.path(out_dir, "rename_key.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

## If outputting full sample names, we just need to subset the input file
if (out_names == "full") {
  system(sprintf("bcftools view %s -S %s -O %s -o %s",
                 in_vcf_file,
                 file.path(out_dir, "samps_subset.txt"),
                 out_fmt,
                 out_vcf_file))

## Otherwise we must first subset, and then reformat the sample names
} else {
  system(sprintf("bcftools view %s -S %s -O %s -o %s",
                 in_vcf_file,
                 file.path(out_dir, "samps_subset.txt"),
                 out_fmt,
                 tmp_var_file))
  system(sprintf("bcftools reheader %s -s %s -o %s",
                 tmp_var_file,
                 file.path(out_dir, "rename_key.txt"),
                 out_vcf_file))
  file.remove(tmp_var_file)
}

## File cleanup
file.remove(file.path(out_dir, "samps_subset.txt"))
file.remove(file.path(out_dir, "rename_key.txt"))
file.remove(file.path(out_dir, "PSC.txt"))
