# Generate DEGgo heatmap

Generates a clustered heatmap of top differentially expressed genes.

## Usage

``` r
plot_heatmap(
  vsd,
  res_df,
  metadata,
  contrast = NULL,
  sample_subset = NULL,
  metadata_filter = NULL,
  top_n_heatmap = 50,
  padj_cutoff = 0.05,
  main = "Top Differentially Expressed Genes",
  output_dir = "DEGgo_out",
  filename = "Heatmap",
  fallback = TRUE,
  annotation_cols = c("condition", "treatment", "sex", "tissue"),
  annotation_colors = NULL,
  order_by = NULL,
  scale_rows = TRUE,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 12,
  fontsize_col = 10,
  width = 8,
  height = 10
)
```

## Arguments

- vsd:

  Variance-stabilized dataset.

- res_df:

  Differential expression results data frame.

- metadata:

  Sample metadata data frame.

- contrast:

  Contrast name to plot.

- sample_subset:

  Optional sample vector to retain.

- metadata_filter:

  Optional named list used to filter metadata.

- top_n_heatmap:

  Number of genes displayed in the heatmap.

- padj_cutoff:

  Adjusted p-value threshold.

- main:

  Heatmap title.

- output_dir:

  Directory for exported plots.

- filename:

  Output file name without extension.

- fallback:

  Logical. If TRUE, use top ranked genes when no significant genes pass
  `padj_cutoff`.

- annotation_cols:

  Metadata columns shown as heatmap annotations.

- annotation_colors:

  Named list of annotation colors.

- order_by:

  Metadata columns used to order samples.

- scale_rows:

  Logical; scale rows before plotting.

- cluster_rows:

  Logical; cluster genes.

- cluster_cols:

  Logical; cluster samples.

- fontsize_row:

  Row label font size.

- fontsize_col:

  Column label font size.

- width:

  Plot width in inches.

- height:

  Plot height in inches.

## Value

Heatmap expression matrix, invisibly.
