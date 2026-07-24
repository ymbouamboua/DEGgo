# Download the GSE98965 baboon transcriptomic atlas

Download the processed gene-level FPKM expression matrix from the Gene
Expression Omnibus accession GSE98965.

## Usage

``` r
download_gse98965(data_dir = "GSE98965")
```

## Arguments

- data_dir:

  Character string. Directory in which the GSE98965 supplementary files
  are stored. The directory is created when it does not exist.

## Value

A `data.frame` containing the complete GSE98965 processed expression
table.

## Details

The dataset corresponds to the baboon diurnal transcriptome atlas
described by Mure et al. and contains gene-level expression values
across multiple tissues and Zeitgeber Time points.

If the expected supplementary file is already present in `data_dir`, it
is reused and is not downloaded again.

This function requires the Bioconductor package GEOquery. Supplementary
files are retrieved with
[`GEOquery::getGEOSuppFiles()`](http://seandavi.github.io/GEOquery/reference/getGEOSuppFiles.md).

## References

Mure LS et al. Diurnal transcriptome atlas of a primate across major
neural and peripheral tissues. Science. 2018.

GEO accession: GSE98965.

## Examples

``` r
if (FALSE) { # \dontrun{
baboon_expression <- download_gse98965(
  data_dir = "GSE98965"
)

dim(baboon_expression)
} # }
```
