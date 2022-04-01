#!/bin/bash

## Export multi-sample VCF/BCF from TileDB
##
## This script will export a multi-sample VCF or BCF file from a TileDB instance.
## As always, writing will be faster for a BCF file, though a corresponding
## VCF file may take up less disk space.
##
## The script takes as input the path to the TileDB instance, the path to the
## output VCF/BCF file, the path to a text file listing samples to export (one
## sample per line), and the chunk size, or
## maximum number of samples to export at once. Note that the chunk size should
## not be greater than the system-imposed maximum number of files that may be
## simultaneously opened. Check this on Linux with ulimit -n. Usually it is set
## to 1,024.
##
## The user can optionally specify a .bed file to only export certain genomic
## regions. Set bed_file to some string that isn't a valid file name to disable 
## this feature.
################################################################################


#### User-Defined Constants ###########

db_path="/autofs/bioinformatics-ward/tiledb_test/500samps"
vcf_file="/autofs/bioinformatics-ward/tiledb_test/exported/all_500_merged.bcf"
samps_file="/autofs/bioinformatics-ward/tiledb_test/500samp_names.txt"
bed_file="none"
chunk_size=100


#### Executable #######################

## Get output format from vcf_file extension
ext="${vcf_file##*.}"
if [[ $ext == "gz" ]]; then
    out_fmt="z"
elif [[ $ext == "bcf" ]]; then
    out_fmt="b"
else
    echo
    echo "Please supply either a .vcf.gz or .bcf file for vcf_file"
    exit 1;
fi

## Create a temp. directory in the directory of the output VCF/BCF
out_dir=$(dirname "$vcf_file")
cvcf_tmp=$(mktemp -d -p "$out_dir")
ssvcf_tmp=$(mktemp -d -p "$out_dir")


## Pretty simple if we request few enough samples to be processed all at once
## just export, index, and merge
n_samps=$(wc -l < "$samps_file")
if [[ $n_samps -lt $chunk_size ]]; then
    if [[ -f "$bed_file" ]]; then
        tiledbvcf export --uri "$db_path" \
            --output-format b \
            --samples-file "$samps_file" \
            --regions-file "$bed_file" \
            --log-level warn \
            --output-dir "$ssvcf_tmp"
    else
        tiledbvcf export --uri "$db_path" \
            --output-format b \
            --samples-file "$samps_file" \
            --log-level warn \
            --output-dir "$ssvcf_tmp"
    fi
    for f in "$ssvcf_tmp"/*.bcf; do bcftools index -c "$f"; done
    bcftools merge --no-version "$ssvcf_tmp"/*.bcf -Ou |
        bcftools view - -m 2 -O "$out_fmt" -o "$vcf_file"
    
## It's more involved if we request more samples than the number of files
## the system can have open at once. In this case we need to split in two stages
else
    ## Figure out total number of chunks we need to split samples into
    ## Note that bash only performs integer division - remainder (if any) is removed
    max_ind=$(($n_samps / $chunk_size))
    n_chunks=$(($max_ind + 1))

    ## For each chunk, output individual sample BCFs and merge together
    for i in $(seq 0 $max_ind); do
        if [[ -f "$bed_file" ]]; then
            tiledbvcf export --uri "$db_path" \
                --output-format b \
                --samples-file "$samps_file" \
                --regions-file "$bed_file" \
                --log-level warn \
                --sample-partition "${i}:${n_chunks}" \
                --output-dir "$ssvcf_tmp"
        else
            tiledbvcf export --uri "$db_path" \
                --output-format b \
                --samples-file "$samps_file" \
                --log-level warn \
                --sample-partition "${i}:${n_chunks}" \
                --output-dir "$ssvcf_tmp"
        fi
        for f in "$ssvcf_tmp"/*.bcf; do bcftools index -c "$f"; done
        bcftools merge --no-version "$ssvcf_tmp"/*.bcf -Ob -o "${cvcf_tmp}/${i}.bcf"
        rm "$ssvcf_tmp"/*
    done

    ## Perform final, secondary merge
    for f in "$cvcf_tmp"/*.bcf; do bcftools index -c "$f"; done
    bcftools merge --no-version "$cvcf_tmp"/*.bcf -Ou |
        bcftools view - -m 2 -O "$out_fmt" -o "$vcf_file"
fi

bcftools index -c "$vcf_file"

rm -rf "$cvcf_tmp"
echo "Removed temp directory $cvcf_tmp"
rm -rf "$ssvcf_tmp"
echo "Removed temp directory $ssvcf_tmp"
