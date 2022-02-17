#!/bin/bash

## Ingest (import) VCF/BCF file into TileDB
##


#### User-Defined Constants ###########

vcf_file="/autofs/bioinformatics-ward/2021_Norgrains_IL_merged_VCF/filt_80miss_3maf_10het/all_regions_samp_filt.bcf"
db_path="/autofs/bioinformatics-ward/norgrains_gbs_tiledb"
chunk_size=500


#### Executable #######################

## Create a temp. directory in the DB's parent directory
db_dir=$(dirname "$db_path")
samps_tmp=$(mktemp -d -p "$db_dir")
cvcf_tmp=$(mktemp -d -p "$db_dir")
ssvcf_tmp=$(mktemp -d -p "$db_dir")

## Get list of samples from VCF/BCF
bcftools query -l "$vcf_file" > "${samps_tmp}/samps.txt"
n_samps=$(wc -l < "${samps_tmp}/samps.txt")

## Pretty simple if the VCF has few enough samples to be processed all at once
## just split into individual samples, index, and import
if [[ $n_samps -lt $chunk_size ]]; then
    bcftools +split "$vcf_file" -Ob -o "$ssvcf_tmp"
    for f in "$ssvcf_tmp"/*.bcf; do bcftools index -c "$f"; done
    tiledbvcf store --uri "$db_path" --log-level warn "$ssvcf_tmp"/*.bcf

## It's much more involved if the VCF has more samples than the number of files
## the system can have open at once. In this case we need to split in two stages
else

    ## Split up our samples file into chunks. For each we create a 3-col file
    ## with sample name, new sample name (just dash to leave as-is) and the
    ## prefix of the BCF file to output to
    split -l $chunk_size "${samps_tmp}/samps.txt" "${samps_tmp}/"
    rm "${samps_tmp}/samps.txt"
    for f in "${samps_tmp}"/*; do
        base=$(basename "$f")
        awk -v fpref="$base" 'BEGIN {OFS = "\t"}{print $0, "-", fpref}' "$f" > "${f}_update.txt"
    done

    ## Concatenate back together and perform first split
    cat "${samps_tmp}"/*update.txt > "${samps_tmp}/split_key.txt"
    bcftools +split "$vcf_file" -G "${samps_tmp}/split_key.txt" -Ob -o "$cvcf_tmp"

    ## Now loop through each primary split BCF and perform the secondary split
    ## (by individual samples); upload all to the DB
    for f in "$cvcf_tmp"/*.bcf; do
        bcftools +split "$f" -Ob -o "$ssvcf_tmp"
        for g in "$ssvcf_tmp"/*.bcf; do bcftools index -c "$g"; done
        tiledbvcf store --uri "$db_path" --log-level warn "$ssvcf_tmp"/*.bcf
        rm "$ssvcf_tmp"/*
    done
fi


rm -rf "$samps_tmp"
echo "Removed temp directory $samps_tmp"
rm -rf "$cvcf_tmp"
echo "Removed temp directory $cvcf_tmp"
rm -rf "$ssvcf_tmp"
echo "Removed temp directory $ssvcf_tmp"
