# Run sample-level RNA-seq quality control

Generate sample-level quality control outputs from a raw or filtered
count matrix, including a Spearman sample correlation heatmap,
hierarchical clustering dendrogram, and sample correlation matrix.

## Usage

``` r
run_sample_qc(
  counts,
  metadata,
  output_dir,
  annotation_cols = c("condition", "treatment", "sex", "tissue"),
  dpi = 300
)
```

## Arguments

- counts:

  Numeric count matrix with genes in rows and samples in columns.

- metadata:

  Sample metadata data frame. Row names must match sample names in
  `colnames(counts)`.

- output_dir:

  Output directory where QC files will be saved.

- annotation_cols:

  Character vector of metadata columns used as heatmap annotations.

- dpi:

  Resolution used for the hierarchical clustering PNG.

## Value

Invisibly returns a list containing the sample correlation matrix,
output directory, and paths to generated files.

## Details

This function is designed to be used inside the DEGgo workflow after
count filtering and metadata matching. It saves all QC files directly
into the provided output directory.

The function uses `log2(count + 1)` transformed counts to compute
Spearman sample correlations and hierarchical clustering based on
Euclidean distance.

## Examples

``` r
if (FALSE) { # \dontrun{
qc <- run_sample_qc(
  counts = counts,
  metadata = metadata,
  output_dir = "pairwise_QC"
)
} # }
```
