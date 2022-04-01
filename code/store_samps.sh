#!/bin/bash

## Ingest (import) VCF/BCF file into TileDB
##
## This script will import a multi-sample VCF/BCF file into a TileDB instance.
## As always, processing will be faster using a BCF file, though a corresponding
## VCF file may take up less disk space.
##
## The script takes as input the path to the input VCF/BCF file, the path to 
## the TileDB instance (must already be initialized), and the chunk size, or
## maximum number of samples to import at once. Note that the chunk size should
## not be greater than the system-imposed maximum number of files that may be
## simultaneously opened. Check this on Linux with ulimit -n. Usually it is set
## to 1,024.
##
## The user can optionally specify a .bed file to only import certain genomic
## regions. This isn't necessarily recommended and it's often best to import
## everything and then restrict to specific regions when exporting from the TileDB.
## Set bed_file to some string that isn't a valid file name to disable this feature.
##
## Likewise, the user can optionally specify a text file listing samples (one per
## line). If this file is specified, then only these samples will be subset from
## the input VCF/BCF and added to the DB. Once again, set this to some string
## that isn't a valid file name to disable.
################################################################################


#### User-Defined Constants ###########

vcf_file="/autofs/bioinformatics-ward/2021_Norgrains_IL_merged_VCF/filt_80miss_3maf_10het/all_regions_samp_filt.bcf"
db_path="/autofs/bioinformatics-ward/norgrains_gbs_tiledb"
ref_file=""
samps_file="/autofs/bioinformatics-ward/2021_Norgrains_IL_merged_VCF/filt_80miss_3maf_10het/genos_diff.txt"
bed_file="none"
chunk_size=500


#### Executable #######################

## Create a temp. directory in the DB's parent directory
db_dir=$(dirname "$db_path")
samps_tmp=$(mktemp -d -p "$db_dir")
cvcf_tmp=$(mktemp -d -p "$db_dir")
ssvcf_tmp=$(mktemp -d -p "$db_dir")

## Get list of samples from VCF/BCF
bcftools query -l "$vcf_file" | sort > "${samps_tmp}/samps.txt"
if [[ -f "$samps_file" ]]; then 
    cat "$samps_file" | sort > "${samps_tmp}/usr_samps.txt"
    comm -12 "${samps_tmp}/samps.txt" "${samps_tmp}/usr_samps.txt" > "${samps_tmp}/tmp.txt"
    mv "${samps_tmp}/tmp.txt" "${samps_tmp}/samps.txt"
    rm "${samps_tmp}/usr_samps.txt"
fi
n_samps=$(wc -l < "${samps_tmp}/samps.txt")

## Pretty simple if the VCF has few enough samples to be processed all at once
## just split into individual samples, index, and import
if [[ $n_samps -lt $chunk_size ]]; then
    
    if [[ -f "$bed_file" ]]; then
        bcftools view "$vcf_file" -S "${samps_tmp}/samps.txt" -R "$bed_file" --force-samples -Ou |
            bcftools norm - -f "$ref_file" -c s -m +both -Ou |
            bcftools +split - -Ob -o "$ssvcf_tmp"
    else
        bcftools view "$vcf_file" -S "${samps_tmp}/samps.txt" --force-samples -Ou |
            bcftools norm - -f "$ref_file" -c s -m +both -Ou |
            bcftools +split - -Ob -o "$ssvcf_tmp"
    fi
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
    cut -f 1 "${samps_tmp}/split_key.txt" > "${samps_tmp}/samps.txt"

    if [[ -f "$bed_file" ]]; then
        bcftools view "$vcf_file" -S "${samps_tmp}/samps.txt" -R "$bed_file" --force-samples -Ou |
            bcftools norm - -f "$ref_file" -c s -m +both -Ou |
            bcftools +split - -G "${samps_tmp}/split_key.txt" -Ob -o "$cvcf_tmp"
    else
        bcftools view "$vcf_file" -S "${samps_tmp}/samps.txt" --force-samples -Ou |
            bcftools norm - -f "$ref_file" -c s -m +both -Ou |
            bcftools +split - -G "${samps_tmp}/split_key.txt" -Ob -o "$cvcf_tmp"
    fi
    
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
