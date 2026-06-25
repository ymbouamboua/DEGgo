# Generate volcano plot for differential expression results

Creates a customizable volcano plot for DEGgo differential expression
results. Genes are automatically classified as Up, Down, or Ns based on
log2 fold-change and p-value thresholds.

## Usage

``` r
plot_volcano(
  res_df,
  gene_col = NULL,
  lfc_col = "log2FoldChange",
  pval_col = "pvalue",
  padj_col = "padj",
  use_padj = TRUE,
  logfc_cutoff = 0.25,
  padj_cutoff = 0.05,
  top_n_labels = 40,
  genes_highlight = NULL,
  colors = c(Up = "#740001", Down = "#6497b1", Ns = "gray70"),
  point_size = 1.5,
  alpha = 0.65,
  title = "Volcano Plot",
  legend = TRUE,
  output_dir = "DEGgo_out",
  filename = "Volcano_Plot",
  save_pdf = TRUE,
  width = 10,
  height = 8,
  dpi = 300,
  style = c("classic", "minimal", "bw", "test", "void", "dirty", "gray"),
  txtsize = 12,
  label_size = 3,
  highlight_label_size = 3.5,
  label_bg = "grey95",
  label_bg_r = 0.12,
  label_force = 1.2,
  label_force_pull = 0.4,
  label_max_overlaps = Inf,
  seed = 42
)
```

## Arguments

- res_df:

  Differential expression result data frame.

- gene_col:

  Column containing gene labels. If unavailable, the function
  automatically tries `"SYMBOL"`, `"gene"`, then `"ENSEMBL"`.

- lfc_col:

  Column containing log2 fold-change values.

- pval_col:

  Column containing raw p-values.

- padj_col:

  Column containing adjusted p-values.

- use_padj:

  Logical. If `TRUE`, use `padj_col`; otherwise use `pval_col`.

- logfc_cutoff:

  Absolute log2 fold-change threshold.

- padj_cutoff:

  P-value or adjusted p-value threshold.

- top_n_labels:

  Number of top significant genes to label.

- genes_highlight:

  Optional character vector of genes to highlight.

- colors:

  Named vector of colors for `Up`, `Down`, and `Ns`.

- point_size:

  Point size.

- alpha:

  Point transparency.

- title:

  Plot title.

- legend:

  Logical. Show legend.

- output_dir:

  Output directory. If not `NULL`, the plot is saved.

- filename:

  Output file name without extension.

- save_pdf:

  Logical. Also save a PDF version.

- width:

  Plot width.

- height:

  Plot height.

- dpi:

  PNG resolution.

- style:

  Plot style passed to `.deggo_theme()`.

- txtsize:

  Base text size.

- label_size:

  Size of standard gene labels.

- highlight_label_size:

  Size of highlighted gene labels.

- label_bg:

  Fill color for label backgrounds.

- label_bg_r:

  Radius of label background boxes.

- label_force:

  Force parameter passed to ggrepel.

- label_force_pull:

  Pull force parameter passed to ggrepel.

- label_max_overlaps:

  Maximum number of overlapping labels.

- seed:

  Random seed for reproducible label placement.

## Value

A `ggplot` object.
