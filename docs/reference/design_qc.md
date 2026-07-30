# Experimental design quality control for bulk RNA-seq

Evaluates metadata structure, repeated experimental units,
PCA-associated effects, optional PERMANOVA, and recommends a
differential-expression model.

## Usage

``` r
design_qc(
  counts,
  metadata,
  sample_col,
  group_col,
  unit_col = NULL,
  batch_cols = NULL,
  gene_col = NULL,
  strip_ensembl_version = TRUE,
  aggregate_duplicates = TRUE,
  transform = c("vst", "rlog", "logcpm"),
  n_top = 500L,
  n_pcs = 5L,
  run_permanova = TRUE,
  permutations = 999L,
  dominance_ratio = 2,
  min_effect_r2 = 0.1,
  verbose = TRUE
)
```

## Arguments

- counts:

  Integer-like matrix/data.frame with genes in rows and samples in
  columns.

- metadata:

  Data frame with one row per sample.

- sample_col:

  Metadata column containing sample identifiers.

- group_col:

  Main biological condition to test.

- unit_col:

  Independent experimental unit (e.g. cage, donor, patient). If NULL, a
  conservative name-based detector is used.

- batch_cols:

  Optional technical/biological covariates to assess.

- gene_col:

  Optional character string specifying the column in `counts` that
  contains gene identifiers. When `NULL`, gene identifiers are expected
  to be stored in the row names.

- strip_ensembl_version:

  Logical. If `TRUE`, remove Ensembl version suffixes from gene
  identifiers, for example converting `"ENSG00000141510.18"` to
  `"ENSG00000141510"`.

- aggregate_duplicates:

  Logical. If `TRUE`, rows with duplicated gene identifiers after
  preprocessing are aggregated by summing their counts. If `FALSE`,
  duplicated identifiers should trigger an error or be retained
  according to the function implementation.

- transform:

  PCA transformation: "vst", "rlog", or "logcpm".

- n_top:

  Number of most variable genes used for PCA.

- n_pcs:

  Number of PCs used for structure assessment.

- run_permanova:

  Run vegan::adonis2 when vegan is installed.

- permutations:

  Number of PERMANOVA permutations.

- dominance_ratio:

  Minimum weighted PCA R2 ratio for declaring one factor stronger than
  the main group.

- min_effect_r2:

  Minimum weighted PCA R2 for an effect to be considered relevant.

- verbose:

  Print a concise summary.

## Value

An object of class "deggo_design_qc".

## Examples

``` r
if (FALSE) { # \dontrun{
qc <- design_qc(
  counts = counts,
  metadata = metadata,
  sample_col = "sample",
  group_col = "Group",
  unit_col = "Cage"
)
print(qc)
plot(qc)
} # }
```
