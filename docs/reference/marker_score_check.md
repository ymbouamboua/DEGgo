# Validate sample identity using marker gene signatures

Computes signature scores from predefined marker gene sets, predicts the
most likely sample identity, and flags potential sample swaps or
annotation mismatches.

## Usage

``` r
marker_score_check(
  counts,
  metadata,
  marker_sets,
  sample_col = "sample",
  group_col = "tissue",
  feature_col = "gene_name",
  log_transform = TRUE,
  plot = TRUE,
  style = "'classic",
  txtsize = 12
)
```

## Arguments

- counts:

  A data frame or matrix containing gene expression counts. Rows
  correspond to genes and columns correspond to samples.

- metadata:

  A data frame containing sample metadata.

- marker_sets:

  Named list of marker gene vectors. Each list element defines a
  biological group or tissue signature.

- sample_col:

  Column name in `metadata` containing sample IDs. Default is
  `"sample"`.

- group_col:

  Column name in `metadata` containing the expected sample annotation
  (e.g. tissue, condition). Default is `"tissue"`.

- feature_col:

  Column name in `counts` containing gene symbols. Default is
  `"gene_name"`.

- log_transform:

  Logical; if `TRUE`, computes scores from `log2(count + 1)` transformed
  values. Default is `TRUE`.

- plot:

  Logical; if `TRUE`, returns a barplot of marker scores. Default is
  `TRUE`.

- style:

  Theme style passed to `.deggo_theme()`.

- txtsize:

  Texte size passed to `.deggo_theme()`

## Value

A list containing:

- scores:

  Data frame with marker scores, predicted groups, and swap flags.

- long:

  Long-format version of the score table.

- swaps:

  Subset of samples where predicted and expected groups differ.

- plot:

  A ggplot object showing marker scores per sample if `plot = TRUE`;
  otherwise `NULL`.

## Details

For each marker set, the function calculates the mean expression of all
detected marker genes across samples. The predicted group corresponds to
the marker set with the highest score.

The function performs a simple signature-based classification:

1.  Matches marker genes against the feature column.

2.  Computes the mean expression of matched genes for each sample.

3.  Predicts the group with the highest signature score.

4.  Flags samples where predicted and expected annotations differ.

Gene matching is case-insensitive.

## Examples

``` r
if (FALSE) { # \dontrun{
marker_sets <- list(
  BAT = c("Ucp1", "Cidea", "Ppargc1a"),
  WAT = c("Adipoq", "Lep", "Fabp4"),
  TESTIS = c("Amh", "Sox9", "Ddx4"),
  OVARY = c("Foxl2", "Fshr", "Bmp15")
)

res <- marker_score_check(
  counts = counts,
  metadata = metadata,
  marker_sets = marker_sets,
  sample_col = "sample",
  group_col = "tissue"
)

head(res$scores)
res$swaps
print(res$plot)
} # }
```
