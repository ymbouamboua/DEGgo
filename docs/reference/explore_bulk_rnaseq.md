# Explore and QC bulk RNA-seq count data

Generates pre-differential-expression quality control plots and metrics
for bulk RNA-seq count data, including library size, detected genes,
expression distribution, density curves, sample correlation,
hierarchical clustering, PCA plots, top variable gene heatmap, optional
marker heatmap, and automatic sample QC flags.

## Usage

``` r
explore_bulk_rnaseq(
  counts,
  metadata,
  gene_col = c("gene_id", "GeneID", "gene", "Gene", "ENSEMBL", "ensembl", "ensembl_id"),
  feature_col = c("gene_name", "SYMBOL", "symbol", "gene_symbol", "external_gene_name"),
  sample_col = c("sample", "Sample", "SAMPLE"),
  annotation_cols = c("condition", "treatment", "sex", "tissue"),
  color_by = NULL,
  markers = c("Adipoq", "Lep", "Ucp1", "Amh", "Amhr2", "Ar", "Akr1c3"),
  top_variable_genes = 500,
  output_dir = "DEGgo_QC",
  width = 8,
  height = 6,
  dpi = 300,
  style = "classic"
)
```

## Arguments

- counts:

  Count table or matrix. Genes should be rows and samples columns. If
  `counts` is a data frame, one column must contain gene identifiers.

- metadata:

  Sample metadata data frame.

- gene_col:

  Candidate gene identifier columns in `counts`.

- feature_col:

  Candidate gene symbol/name columns in `counts`.

- sample_col:

  Candidate sample identifier columns in `metadata`.

- annotation_cols:

  Metadata columns shown as annotations in heatmaps.

- color_by:

  Metadata column used to color library size and density plots. If
  `NULL`, DEGgo tries tissue, condition, treatment, then sex.

- markers:

  Optional marker genes for marker heatmap.

- top_variable_genes:

  Number of top variable genes for heatmap.

- output_dir:

  Output directory.

- width:

  Plot width in inches.

- height:

  Plot height in inches.

- dpi:

  Plot resolution.

- style:

  Theme style passed to `.deggo_theme()`.

## Value

Invisibly returns a list containing QC table, plots, correlation matrix,
PCA plots, top variable genes, marker heatmap object, and output path.

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
