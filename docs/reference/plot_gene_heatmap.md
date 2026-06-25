# Plot expression heatmap for selected genes

Generates a clustered heatmap for user-defined genes.

## Usage

``` r
plot_gene_heatmap(
  counts,
  metadata,
  genes,
  gene_col = c("gene_id", "GeneID", "gene", "Gene", "ENSEMBL", "ensembl", "ensembl_id"),
  feature_col = c("gene_name", "SYMBOL", "symbol", "gene_symbol", "external_gene_name"),
  sample_col = c("sample", "Sample", "SAMPLE"),
  assay_transform = c("log2", "log2cpm"),
  annotation_cols = c("condition", "treatment", "sex", "tissue"),
  annotation_colors = NULL,
  order_by = NULL,
  output_dir = "DEGgo_out",
  filename = "Gene_Expression_Heatmap",
  main = "Selected gene expression heatmap",
  color = (grDevices::colorRampPalette(c("#6497b1", "#F7F7F7", "#740001")))(100),
  breaks = seq(-2, 2, length.out = 101),
  scale_rows = TRUE,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 12,
  fontsize_col = 10,
  width = 10,
  height = 7
)
```

## Arguments

- counts:

  Count matrix or count table.

- metadata:

  Sample metadata.

- genes:

  Genes to display.

- gene_col:

  Gene ID column.

- feature_col:

  Gene symbol column.

- sample_col:

  Sample column.

- assay_transform:

  Expression transformation.

- annotation_cols:

  Metadata columns shown above heatmap.

- annotation_colors:

  Named list of annotation colors.

- order_by:

  Metadata columns used to order samples.

- output_dir:

  Output directory.

- filename:

  Output filename.

- main:

  Heatmap title.

- color:

  Heatmap color palette.

- breaks:

  Numeric vector of color breaks.

- scale_rows:

  Scale genes.

- cluster_rows:

  Cluster genes.

- cluster_cols:

  Cluster samples.

- fontsize_row:

  Row names size.

- fontsize_col:

  Column names size.

- width:

  Plot width.

- height:

  Plot height.

## Value

Expression matrix used for plotting.
