# Run DEGgo bulk RNA-seq downstream analysis

Run a complete and automated bulk RNA-seq downstream analysis workflow.

## Usage

``` r
run_deggo(
  counts,
  metadata,
  gene_col = c("gene_id", "GeneID", "gene", "Gene", "ENSEMBL", "ensembl", "ensembl_id"),
  feature_col = c("gene_name", "SYMBOL", "symbol", "gene_symbol", "external_gene_name"),
  sample_col = c("sample", "Sample", "SAMPLE"),
  prepare_input = TRUE,
  raw_qc = TRUE,
  remove_flagged = FALSE,
  qc_markers = NULL,
  marker_sets = NULL,
  qc_sample_col = NULL,
  qc_output_prefix = "DEGgo_QC",
  output_dir = NULL,
  padj_cutoff = 0.05,
  logfc_cutoff = 0.25,
  top_n_heatmap = 50,
  top_n_labels = 10,
  ontology = c("BP", "MF", "CC"),
  organism = c("human", "mouse", "rat", "custom"),
  orgdb = NULL,
  method = c("DESeq2", "edgeR", "limma"),
  analysis_mode = c("single", "pairwise"),
  contrast = NULL,
  design_formula = ~condition,
  pairwise_group_cols = NULL,
  pairwise_contrast_col = "comparison_group",
  pairwise_contrasts = NULL,
  filter_method = c("count", "cpm", "none"),
  pairwise_mode = c("all", "within_first", "within_second"),
  min_count = 5,
  min_samples = 2,
  min_total = 10,
  generate_report = TRUE,
  report_formats = "html",
  report_template = NULL,
  generate_pptx = FALSE,
  pptx_file = NULL,
  save_reproducibility = TRUE,
  save_clean_inputs = TRUE,
  txtsize = 12,
  seed = 123
)
```

## Arguments

- counts:

  Raw count table, data frame or matrix. Rows should represent
  genes/features and columns should represent samples. If
  `prepare_input = TRUE`, DEGgo attempts to detect gene and sample
  columns automatically using `gene_col`, `feature_col`, and
  `sample_col`.

- metadata:

  Sample metadata data frame. Must contain a sample identifier column
  matching the count matrix sample names. For single analysis, it must
  contain `condition` unless this is encoded in the design/contrast
  workflow. For pairwise analysis, `condition` can be automatically
  generated from `pairwise_group_cols`.

- gene_col:

  Character vector of possible gene identifier column names. Used during
  input preparation.

- feature_col:

  Character vector of possible feature/gene symbol column names. Used
  during input preparation and marker-based checks.

- sample_col:

  Character vector of possible sample identifier column names in
  `metadata`.

- prepare_input:

  Logical. If `TRUE`, automatically prepares and matches the count table
  and metadata using `prepare_counts_metadata()`. If `FALSE`, `counts`
  must already be a numeric matrix with sample names matching metadata
  row names or a sample column.

- raw_qc:

  Logical. If `TRUE`, run exploratory raw sample QC before input
  preparation and differential expression analysis using
  [`explore_bulk_rnaseq()`](https://ymbouamboua.github.io/DEGgo/reference/explore_bulk_rnaseq.md).

- remove_flagged:

  Logical. If `TRUE`, remove samples flagged by the raw QC table before
  continuing the workflow. Default is `FALSE` to keep sample exclusion
  under user control.

- qc_markers:

  Optional character vector of marker genes used during raw and clean
  QC.

- marker_sets:

  Optional named list of marker gene sets used for marker-based
  biological validation with
  [`marker_score_check()`](https://ymbouamboua.github.io/DEGgo/reference/marker_score_check.md).

- qc_sample_col:

  Optional sample column used by
  [`remove_flagged_samples()`](https://ymbouamboua.github.io/DEGgo/reference/remove_flagged_samples.md).
  If `NULL`, the first value of `sample_col` is used.

- qc_output_prefix:

  Character prefix used for QC output directories.

- output_dir:

  Output directory. If `NULL`, results are written to `"DEGgo_results"`.
  A dated DEGgo subdirectory is created automatically.

- padj_cutoff:

  Adjusted p-value cutoff used to define significant DEGs.

- logfc_cutoff:

  Absolute log2 fold-change cutoff used to define significant DEGs.

- top_n_heatmap:

  Number of top genes to display in DEG heatmaps.

- top_n_labels:

  Number of top genes to label in volcano plots.

- ontology:

  Gene Ontology namespace. One of `"BP"`, `"MF"`, or `"CC"`.

- organism:

  Organism used for gene annotation. One of `"human"`, `"mouse"`,
  `"rat"`, or `"custom"`.

- orgdb:

  Optional AnnotationDbi OrgDb object. Required when
  `organism = "custom"`.

- method:

  Differential expression method. One of `"DESeq2"`, `"edgeR"`, or
  `"limma"`. Pairwise mode currently supports DESeq2 only.

- analysis_mode:

  Analysis mode. `"single"` runs one differential expression analysis.
  `"pairwise"` runs multiple pairwise DESeq2 contrasts.

- contrast:

  Optional contrast passed to the differential expression engine in
  single-analysis mode.

- design_formula:

  Design formula used by the differential expression engine, for example
  `~ condition` or `~ batch + condition`.

- pairwise_group_cols:

  Character vector of metadata columns used to define pairwise groups.
  Required when `analysis_mode = "pairwise"`.

- pairwise_contrast_col:

  Name of the metadata column used to store or generate pairwise
  contrast groups.

- pairwise_contrasts:

  Optional named list of pairwise contrasts. If provided, these
  contrasts are used instead of automatically generated contrasts.

- filter_method:

  Low-expression filtering method. One of `"count"`, `"cpm"`, or
  `"none"`.

- pairwise_mode:

  Pairwise comparison mode. One of `"all"`, `"within_first"`, or
  `"within_second"`.

- min_count:

  Minimum count threshold used by count-based filtering.

- min_samples:

  Minimum number of samples required to pass the count/CPM threshold.

- min_total:

  Minimum total count required for a gene to be retained.

- generate_report:

  Logical. If `TRUE`, generate DEGgo report files.

- report_formats:

  Character vector of report formats, for example `"html"`, `"pdf"`, or
  `c("html", "pdf")`.

- report_template:

  Optional path to a custom R Markdown report template. If `NULL`, the
  package DEGgo template is used.

- generate_pptx:

  Logical. If `TRUE`, generate a PowerPoint report.

- pptx_file:

  Optional path to the output PowerPoint file. If `NULL`,
  `"DEGgo_Report.pptx"` is written inside `output_dir`.

- save_reproducibility:

  Logical. If `TRUE`, save reproducibility files, including run objects
  and session information.

- save_clean_inputs:

  Logical. If `TRUE`, save cleaned count and metadata tables after
  sample matching and filtering.

- txtsize:

  Base text size used in DEGgo visualizations.

- seed:

  Random seed for reproducibility.

## Value

A named list containing differential expression results, significant
DEGs, summary tables, plots, GO enrichment results, QC outputs, cleaned
counts and metadata, report paths, PowerPoint path, output directories,
run parameters, manifest, and DEGgo version.

## Details

DEGgo can perform optional raw sample quality control, optional removal
of flagged samples, marker-based biological validation, input
preparation, sample matching, gene identifier cleaning, low-expression
filtering, differential expression analysis, result annotation, PCA,
heatmap and volcano visualization, Gene Ontology enrichment,
reproducibility export, and automated HTML/PDF/PPTX report generation.

The function supports either a single differential expression analysis
or multiple pairwise contrasts. Pairwise mode currently uses DESeq2.

In single mode, DEGgo runs one differential expression analysis using
the selected method and design formula. In pairwise mode, DEGgo builds
or uses predefined pairwise contrasts from metadata columns and runs
DESeq2 for each contrast.

The raw QC step is performed before input preparation and filtering. The
sample QC step is performed after gene filtering and uses the final
matched count matrix and metadata.

Automatic sample removal is disabled by default because exclusion of
samples should usually be reviewed by the analyst.

## See also

`prepare_counts_metadata`, `explore_bulk_rnaseq`,
`remove_flagged_samples`, `marker_score_check`, `run_de`,
`run_deseq2_pairwise`, `run_go_enrichment`, `generate_deggo_report`,
`generate_deggo_pptx`

## Examples

``` r
if (FALSE) { # \dontrun{
results <- run_deggo(
  counts = counts,
  metadata = metadata,
  organism = "mouse",
  method = "DESeq2",
  analysis_mode = "single",
  design_formula = ~ condition,
  contrast = c("condition", "treated", "control"),
  sample_col = "sample"
)

pairwise_results <- run_deggo(
  counts = counts,
  metadata = metadata,
  organism = "mouse",
  method = "DESeq2",
  analysis_mode = "pairwise",
  pairwise_group_cols = c("treatment", "sex", "tissue"),
  pairwise_contrasts = pairwise_contrasts,
  sample_col = "sample",
  raw_qc = TRUE,
  remove_flagged = FALSE
)
} # }
```
