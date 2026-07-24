# Plot a DEGgo Differential Expression Heatmap

Plot a DEGgo Differential Expression Heatmap

## Usage

``` r
plot_heatmap(
  vsd,
  res_df,
  metadata,
  contrast = NULL,
  sample_subset = NULL,
  metadata_filter = NULL,
  top_n_heatmap = 20,
  padj_cutoff = 0.05,
  logfc_cutoff = 0.25,
  main = "Top Differentially Expressed Genes",
  output_dir = "DEGgo_out",
  filename = "Heatmap",
  fallback = TRUE,
  annotation_colors = NULL,
  annotation_cols = NULL,
  order_by = NULL,
  order_levels = NULL,
  scale_rows = TRUE,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = NULL,
  fontsize_col = NULL,
  show_dendrogram = FALSE,
  width = NULL,
  height = NULL,
  show_rownames = NULL,
  show_colnames = NULL,
  palette = "default"
)
```

## Arguments

- vsd:

  Variance-stabilized expression object.

- res_df:

  Differential-expression table, preferably a cleaned table from
  `de_results$sig_deg_clean`.

- metadata:

  Sample metadata.

- contrast:

  Optional contrast vector.

- sample_subset:

  Optional samples to retain.

- metadata_filter:

  Optional named metadata filtering list.

- top_n_heatmap:

  Number of genes displayed.

- padj_cutoff:

  Adjusted P-value cutoff.

- logfc_cutoff:

  Absolute log2 fold-change cutoff used when selecting significant genes
  for the heatmap.

- main:

  Heatmap title.

- output_dir:

  Output directory.

- filename:

  Output filename without extension.

- fallback:

  Use ranked genes if no gene passes the significance cutoff.

- annotation_colors:

  Optional annotation-color list.

- annotation_cols:

  Metadata annotation columns.

- order_by:

  Metadata variables used to order samples.

- order_levels:

  Optional named list defining sample-group order.

- scale_rows:

  Scale expression by gene.

- cluster_rows:

  Cluster genes.

- cluster_cols:

  Cluster samples.

- fontsize_row:

  Row-label size.

- fontsize_col:

  Column-label size.

- show_dendrogram:

  Display dendrograms.

- width:

  Output width.

- height:

  Output height.

- show_rownames:

  Display row names.

- show_colnames:

  Display column names.

- palette:

  DEGgo palette name.

## Value

Invisibly returns the plotted expression matrix.
