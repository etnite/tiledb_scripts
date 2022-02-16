#!/bin/bash

## Export multi-sample VCF/BCF from TileDB
##

#### User-Defined Constants ###########

db_path="/autofs/bioinformatics-ward/tiledb_test/500samps"
vcf_file="/autofs/bioinformatics-ward/tiledb_test/exported/all_500_merged.bcf"
samps_file="/autofs/bioinformatics-ward/tiledb_test/500samp_names.txt"
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
    tiledbvcf export --uri "$db_path" \
        --output-format b \
        --samples-file "$samps_file" \
        --log-level info \
        --output-dir "$ssvcf_tmp"
    for f in "$ssvcf_tmp"/*.bcf; do bcftools index -c "$f"; done
    bcftools merge "$ssvcf_tmp"/*.bcf -O "$out_fmt" -o "$vcf_file"
    
## It's more involved we request more samples than the number of files
## the system can have open at once. In this case we need to split in two stages
else
    ## Figure out total number of chunks we need to split samples into
    ## Note that bash only performs integer division - remainder (if any) is removed
    n_chunks=$(($n_samps / $chunk_size))
    n_chunks=$(($n_chunks + 1))

    ## For each chunk, output individual sample BCFs and merge together
    for i in $(seq 1 $n_chunks); do
        tiledbvcf export --uri "$db_path" \
            --output-format b \
            --samples-file "$samps_file" \
            --log-level info \
            --sample-partition "${i}:${n_chunks}" \
            --output-dir "$ssvcf_tmp"
        for f in "$ssvcf_tmp"/*.bcf; do bcftools index -c "$f"; done
        bcftools merge "$ssvcf_tmp"/*.bcf -Ob -o "${cvcf_tmp}/${i}.bcf"
        rm "$ssvcf_tmp"/*
    done

    ## Perform final, secondary merge
    for f in "$cvcf_tmp"/*.bcf; do bcftools index -c "$f"; done
    bcftools merge "$cvcf_tmp"/*.bcf -O "$out_fmt" -o "$vcf_file"
fi


#rm -rf "$cvcf_tmp"
#echo "Removed temp directory $cvcf_tmp"
#rm -rf "$ssvcf_tmp"
#echo "Removed temp directory $ssvcf_tmp"