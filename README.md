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

TileDB-VCF seems to be geared somewhat more towards a GATK workflow, where
single-sample gVCF files are created prior to SNP calling. Typically I avoid
using GATK because, as of this writing, it still doesn't support the newer .csi
indices required for larger chromosomes (e.g. many of the wheat chromosomes).

Therefore my use case is a little roundabout. We will typically call variants
using BCFTools or Freebayes. Both of these tools create cohort or multi-sample
VCF/BCF files. These must then be split into single-sample VCF/BCF files using
`bcftools +split` prior to ingestion into the DB. On the other end, single sample
VCF/BCF files are exported from the DB, then combined using `bcftools merge`. This
workflow is admittedly a bit cumbersome, though apparently multi-sample import and
export are on the project's roadmap.

### Database size

A typical TileDB for VCF data will be larger than a corresponding single multi-sample
VCF/BCF file; in my testing about 4X or 5X larger, though this heavily depends upon 
the sparseness of the data. This is because the DB stores metadata (i.e. header data) for each sample, whereas in a 
multi-sample file this information is stored only once. On the other hand this
approach allows for much better record keeping regarding how each particular
sample was processed.

A TileDB will be substantially smaller than the sum of the single-sample files
it is built out of. This is because, unlike a flat file, the TileDB does
not store missing data.

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