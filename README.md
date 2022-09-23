# tiledb_scripts

A set of scripts for setting up and using a local TileDB-VCF instance with the CLI

See the official documentation and Git repository at:

https://docs.tiledb.com/main/integrations-and-extensions/population-genomics/

https://github.com/TileDB-Inc/TileDB-VCF

This repository houses code to assist with the two basic tasks when working with
a TileDB-VCF instance:

1. "Ingesting" (i.e. importing) single-sample VCF/BCF files into the DB
2. Exporting single-sample VCF/BCF files out of the DB

While TileDB can, for the time being, only handle import/export of single-sample
files, the scripts in this repository are designed to help facilitate import/export
of multi-sample/cohort files.

## Considerations

While I believe that a system like TileDB is a very sensible alternative to
simply storing genetic variants in giant VCF/BCF files, there are a few pros, cons, and general
considerations to keep in mind:

### Updating

The main impetus for using TileDB is that new samples can be added to the DB
without modifying the existing DB contents. This avoids having to rerun SNP calling
any time new data is added. One thing to note is that, as far as I can tell,
there is no way to remove a sample once it has been added to the DB.

### Front-loaded calling

One advantage of using TileDB-VCF, as opposed to GATK's GenomicsDB, is that the 
variant calling process can happen prior to DB import depending on the software used. 
Therefore, once new samples are imported into the DB, generating a new
multi-sample VCF/BCF file entails an export step followed by a merging step(s),
which are much faster.

### Workflow

My use case is a little roundabout as I don't usually create single-sample VCF
files as part of my variant calling process. I will typically call variants
using BCFTools or Freebayes. Both of these tools create cohort or multi-sample
VCF/BCF files. These must then be split into single-sample VCF/BCF files using
`bcftools +split` prior to ingestion into the DB. On the other end, single sample
VCF/BCF files are exported from the DB, then combined using `bcftools merge`. This
workflow is admittedly a bit cumbersome, though apparently multi-sample import and
export are on the project's roadmap.

### Database size

A typical TileDB for VCF data will be larger than a corresponding single multi-sample
VCF/BCF file; in my testing about 4X or 5X larger, though this heavily depends upon 
the sparseness of the data. This is because the DB stores metadata (e.g. header data) for each sample, whereas in a 
multi-sample file this information is stored only once. On the other hand this
approach allows for much better record keeping regarding how each particular
sample was processed.

A TileDB will be substantially smaller than the sum of the single-sample files
it is built out of. This is because, unlike a flat file, the TileDB does
not store any missing data.

## Installation

The easiest way to install TileDB-VCF is via Conda (note that this repository
focuses mostly on using the CLI, though there is also a Python package available):

```
# Install the CLI and shared library
conda install -c conda-forge -c bioconda -c tiledb libtiledbvcf

# Install the Python package
conda install -c conda-forge -c bioconda -c tiledb tiledbvcf-py

# Run tiledbvcf --version to check if install worked
tiledbvcf --version
```

## Initializing a new DB

Creating a new, empty DB is very straightforward (note that in all TileDB
commands, the name of the DB is given with `--uri`):

```
tiledbvcf create --uri path/to/db
```

## Storing samples

This repository assumes that the user has a multi-sample VCF or BCF file that
they would like to import into a TileDB-VCF instance. The script to perform this
step is `code/store_samps.sh`. A few notes on this process:

1. It is possible to import a sample that is already present within the DB. In 
this case I think(?) TileDB replaces the sample (unfortunately without a warning).
I would like to add some exception handling for this, but for now it's the user's
responsibility to watch out for sample collisions
2. In the case of importing large files, `bcftools +split` will run into issues
with the system-imposed limit on the number of files that can be open simultaneously.
You can check this with `ulimit -n`. Typically on Linux systems this will be 1,024.
If the number of samples to import is greater than a user-specified threshold,
they will be imported in batches.

## Exporting samples

This repository assumes that the user wants to create a multi-sample VCF/BCF file
from an existing TileDB-VCF instance using the script `code_export_samps.sh`. Once again, for large numbers of samples
`bcftools merge` will run into issues with the number of files that can be open
at once. Therefore in these cases a stepwise multi-level merging is performed.
The export step can also take a .bed file specifying which region(s) to export.

## Filtering

The contents of the TileDB will depend upon any filtering performed on the input
VCFs/BCFs prior to import. Typically I lightly filter raw files by SNPwise
missingness, SNPwise het. call proportion, MAF, and sample-wise het. call
proportion prior to importing, just to get rid of SNPs and samples that won't likely
be very useful. Likewise, the multi-sample VCF/BCF
generated from exported samples can be filtered however the user wants. 

## Dealing with Repeated Samples

One thing we need to consider when using this type of database is how to deal with
samples that have been genotyped multiple times. Currently, we store samples in the form:

```
<library_prep_num>.<project>.<state>.<year>.<genotype>
```

So we could have the same genotype appear multiple times, if it was genotyped in different projects, years, etc.
For instance:

```
241587.NIFA_proj.IL.20.IL00-8530
255817.NUWWSN.IL.22.IL00-8530
```

Once these two repeated samples of the same genotype are in the database, they cannot be merged - merging typically
happens during the calling process (if we so choose). Therefore, the question arises, which sample should we use?
This repository has a few scripts to deal with this situation:

### Find repeated samples

The script `code/find_repeat_samps.R` is a simple script to find which genotypes present in a TileDB instance are
repeated. We first generate a list of everything present in the database:

```
tiledbvcf list --uri /path/to/db > sample_list.txt
```

Then we just supply this list as a positional argument to find_repeat_samps.R:

```
Rscript find_repeat_samps.R sample_list.txt
```

This outputs a text file in the same location, in this case labeled `sample_list_repeats_only.txt`

### Classify repeated samples

Say we have a particular genotype that has been genotyped 8 times in our database. Since we can't merge these together,
we need to find a good sample that represents this group. However, how do we know which of these 8 repeats are legit?
There are a few things we need to careful of:

1. Cases where one or more repeats are outliers and are not as closely related to the others as we assume
2. Cases where none of the repeats are very... repeatable (i.e. none of them appear very similar and there is essentially no consensus)

The script `code/classify_repeat_samps.R` can be used for this purpose. It takes a VCF file as input, and 
calculates an identity-by-state (IBS) matrix for all samples. IBS is essentially a simply matching coefficient
that measures how often alleles differ between two samples. A value of 0 would indicate perfect concurrence
between two samples, while 1 would indicate two samples which share no alleles.

Based on the sample naming scheme above, the script finds repeated genotypes and then classifies them into three
categories in an output .csv file:

1. retain
2. outlier
3. no_consens

The threshold for determining whether or not samples are "similar" or not is a value between 0 and 1, with values closer to
0 being more stringent (I often use 0.1). The input VCF should not have a high level of missing data, as this may
lead to inflated IBS coefficients.

The .csv file produced by `classify_repeat_samps.R` can then be used to remove bad replicates from the VCF file.

### Remove repeats from a VCF

Assuming we have already removed all "bad" repeated samples from a VCF file, we can select one representative sample
for each remaining group of repeats. The script `code/remove_repeats_from_vcf.R` can be used to perform this task.
For each set of repeated samples, it will retain the sample with the lowest amount of missing data.
