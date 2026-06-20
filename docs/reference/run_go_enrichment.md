# Run Gene Ontology enrichment analysis

Performs Gene Ontology enrichment analysis separately on upregulated and
downregulated significant differentially expressed genes using Entrez
identifiers.

## Usage

``` r
run_go_enrichment(
  sig_deg,
  comparison = NULL,
  ontology = c("BP", "MF", "CC"),
  orgdb,
  output_dir = "DEGgo_out",
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.2,
  p_adjust_method = "BH",
  logfc_col = "log2FoldChange",
  entrez_col = "ENTREZID",
  min_genes = 10
)
```

## Arguments

- sig_deg:

  Significant DEG data frame containing Entrez IDs and log2 fold-change
  values.

- comparison:

  Character. Comparison name used to label GO enrichment results and
  exported files.

- ontology:

  Character. GO ontology to test. One of `"BP"`, `"MF"` or `"CC"`.

- orgdb:

  OrgDb annotation database used by
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).

- output_dir:

  Character. Output directory.

- pvalue_cutoff:

  Numeric. Raw p-value cutoff.

- qvalue_cutoff:

  Numeric. q-value cutoff.

- p_adjust_method:

  Character. Multiple testing correction method.

- logfc_col:

  Character. Column containing log2 fold-change values.

- entrez_col:

  Character. Column containing Entrez gene identifiers.

- min_genes:

  Integer. Minimum number of genes required to perform enrichment.

## Value

A named list containing:

- `go_results`:

  GO enrichment results table.

- `ontology`:

  Ontology tested.

- `comparison`:

  Comparison identifier.

- `n_up_genes`:

  Number of upregulated genes tested.

- `n_down_genes`:

  Number of downregulated genes tested.
