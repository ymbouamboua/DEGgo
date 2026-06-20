# Validate raw count table

Performs sanity checks on bulk RNA-seq count tables.

## Usage

``` r
check_raw_counts(
  counts,
  gene_col = c("gene_id", "GeneID", "gene", "Gene", "ENSEMBL", "ensembl"),
  feature_col = c("gene_name", "SYMBOL", "symbol"),
  min_library_size = 1e+06,
  min_detected_genes = 5000
)
```

## Arguments

- counts:

  Raw count table.

- gene_col:

  Possible gene ID columns.

- feature_col:

  Possible gene symbol columns.

- min_library_size:

  Minimum acceptable library size.

- min_detected_genes:

  Minimum acceptable detected genes.

## Value

List of QC results.
