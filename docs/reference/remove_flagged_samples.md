# Remove flagged samples from counts and metadata

Removes samples flagged for exclusion in a QC table and returns
synchronized counts and metadata objects.

## Usage

``` r
remove_flagged_samples(
  counts,
  metadata,
  qc_table,
  sample_col = "sample_original",
  remove_col = "recommend_remove",
  gene_cols = c("gene_id", "gene_name"),
  verbose = TRUE
)
```

## Arguments

- counts:

  Count matrix or data frame. Rows are genes and columns are samples.
  Gene annotation columns are preserved.

- metadata:

  Sample metadata.

- qc_table:

  QC table containing sample identifiers and removal flags.

- sample_col:

  Column containing sample names. Default is `"sample_original"`.

- remove_col:

  Logical column indicating samples to remove. Default is
  `"recommend_remove"`.

- gene_cols:

  Character vector of annotation columns to preserve. Default is
  `c("gene_id", "gene_name")`.

- verbose:

  Logical; print summary messages. Default is `TRUE`.

## Value

A list with counts, metadata, and removed sample names.

## Examples

``` r
if (FALSE) { # \dontrun{
qc <- explore_bulk_rnaseq(counts, metadata)

cleaned <- remove_flagged_samples(
  counts = counts,
  metadata = metadata,
  qc_table = qc$qc
)

counts <- cleaned$counts
metadata <- cleaned$metadata
} # }
```
