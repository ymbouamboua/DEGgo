# Plot Expression Heatmap for Selected Genes

Generate a publication-ready heatmap for user-defined genes using raw
counts and sample metadata. The function matches samples, transforms
expression, optionally scales genes by row, adds metadata annotations,
and exports a PNG heatmap.

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
  color = (grDevices::colorRampPalette(c("#6497B1", "#F7F7F7", "#740001")))(100),
  breaks = seq(-2, 2, length.out = 101),
  scale_rows = TRUE,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = NULL,
  fontsize_col = NULL,
  show_dendrogram = FALSE,
  width = NULL,
  height = NULL,
  show_rownames = NULL,
  show_colnames = NULL
)
```

## Arguments

- counts:

  Count matrix or count table. If a data frame is provided, one column
  must contain gene identifiers.

- metadata:

  Sample metadata data frame.

- genes:

  Character vector of gene IDs or gene symbols to display.

- gene_col:

  Candidate gene ID columns.

- feature_col:

  Candidate gene symbol/name columns.

- sample_col:

  Candidate sample identifier columns in `metadata`.

- assay_transform:

  Expression transformation, either `"log2"` or `"log2cpm"`.

- annotation_cols:

  Metadata columns shown above the heatmap.

- annotation_colors:

  Optional annotation colors passed to `pheatmap`.

- order_by:

  Optional metadata columns used to order samples.

- output_dir:

  Output directory.

- filename:

  Output filename without extension.

- main:

  Heatmap title.

- color:

  Heatmap color palette.

- breaks:

  Numeric vector of color breaks.

- scale_rows:

  Logical. If `TRUE`, scale expression by gene.

- cluster_rows:

  Logical. If `TRUE`, cluster genes.

- cluster_cols:

  Logical. If `TRUE`, cluster samples.

- fontsize_row:

  Optional row label font size. If `NULL`, chosen automatically.

- fontsize_col:

  Optional column label font size. If `NULL`, chosen automatically.

- show_dendrogram:

  Logical. If TRUE, display row/column dendrograms when clustering is
  enabled.

- width:

  Optional plot width in inches. If `NULL`, chosen automatically.

- height:

  Optional plot height in inches. If `NULL`, chosen automatically.

- show_rownames:

  Optional logical. If `NULL`, shown automatically for heatmaps with 80
  genes or fewer.

- show_colnames:

  Optional logical. If `NULL`, shown automatically for heatmaps with 60
  samples or fewer.

## Value

Invisibly returns the expression matrix used for plotting.

## Details

By default, genes are clustered but dendrograms are hidden for a clean
report-friendly DEGgo visualization. Samples are not clustered by
default, preserving metadata-defined order.
