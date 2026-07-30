# Run DEGgo bulk RNA-seq downstream analysis

Automated bulk RNA-seq differential expression workflow including QC,
preprocessing, differential expression, annotation, visualization, GO
enrichment, reporting, and reproducibility exports.

## Usage

``` r
run_deggo(
  counts,
  metadata,
  project_name = NULL,
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
  min_expr_count = 20,
  min_expr_samples = 3,
  min_prevalence = 0.6,
  max_sample_fraction = 0.45,
  max_group_sample_fraction = 0.45,
  min_group_mean = 10,
  min_group_median = 20,
  max_group_cv = NULL,
  expr_filter_groups = "auto",
  clean_deg_tables = TRUE,
  ontology = c("BP", "MF", "CC"),
  organism = c("human", "mouse", "rat", "custom"),
  orgdb = NULL,
  method = c("DESeq2", "edgeR", "limma", "dream"),
  dream_ddf = c("adaptive", "Satterthwaite", "Kenward-Roger"),
  dream_n_cores = 1L,
  analysis_mode = c("single", "pairwise"),
  contrast = NULL,
  design_formula = ~condition,
  pairwise_group_cols = NULL,
  pairwise_contrast_col = "group",
  pairwise_contrasts = NULL,
  filter_method = c("count", "cpm", "none"),
  pairwise_mode = c("all", "within_first", "within_second"),
  min_count = 5,
  min_samples = 2,
  min_total = 10,
  rhythmicity_analysis = FALSE,
  rhythmicity_time_col = "time",
  rhythmicity_group_col = NULL,
  rhythmicity_assay = c("vst", "normalized", "log2_normalized", "raw"),
  rhythmicity_methods = c("meta2d", "cosinor"),
  rhythmicity_period_range = c(20, 28),
  rhythmicity_cycle_length = 24,
  rhythmicity_cycMethod = c("ARS", "JTK", "LS"),
  cosinor_engine = c("auto", "package", "manual"),
  rhythmicity_plots = TRUE,
  rhythmicity_n_top_plots = 20,
  generate_report = TRUE,
  report_formats = "html",
  report_template = NULL,
  generate_pptx = FALSE,
  pptx_file = NULL,
  save_reproducibility = TRUE,
  save_clean_inputs = TRUE,
  txtsize = 12,
  heatmap_annotation_cols = "auto",
  palette = "default",
  seed = 123
)
```

## Arguments

- counts:

  Raw count matrix or data frame.

- metadata:

  Sample metadata data frame.

- project_name:

  Optional project name shown in HTML, PDF and PowerPoint reports.

- gene_col:

  Candidate gene identifier column names.

- feature_col:

  Candidate gene symbol/name column names.

- sample_col:

  Candidate sample identifier column names.

- prepare_input:

  Logical. Prepare and match input tables.

- raw_qc:

  Logical. Run raw count QC.

- remove_flagged:

  Logical. Remove QC-flagged samples.

- qc_markers:

  Optional marker genes for QC.

- marker_sets:

  Optional named marker gene sets.

- qc_sample_col:

  Optional sample column for QC.

- qc_output_prefix:

  QC output prefix.

- output_dir:

  Output directory.

- padj_cutoff:

  Adjusted p-value cutoff.

- logfc_cutoff:

  Absolute log2 fold-change cutoff.

- top_n_heatmap:

  Number of genes shown in heatmaps.

- top_n_labels:

  Number of genes labelled in volcano plots.

- min_expr_count:

  Minimum expression count for clean DEG filtering.

- min_expr_samples:

  Minimum samples passing expression threshold.

- min_prevalence:

  Minimum prevalence for clean DEG filtering.

- max_sample_fraction:

  Maximum single-sample fraction.

- max_group_sample_fraction:

  Maximum fraction of group expression contributed by a single sample.

- min_group_mean:

  Minimum group mean expression.

- min_group_median:

  Optional minimum median expression required in at least one comparison
  group during post-DE filtering.

- max_group_cv:

  Optional maximum within-group coefficient of variation used during
  post-DE filtering.

- expr_filter_groups:

  Grouping variables for expression filtering.

- clean_deg_tables:

  Logical. Apply post-DEG expression cleaning.

- ontology:

  GO ontology.

- organism:

  Organism name.

- orgdb:

  Optional custom OrgDb object.

- method:

  Differential expression method.

- dream_ddf:

  Méthode utilisée pour calculer les degrés de liberté dans les modèles
  `dream`. L'une de `"adaptive"`, `"Satterthwaite"` ou
  `"Kenward-Roger"`.

- dream_n_cores:

  Nombre de cœurs utilisés par `dream`.

- analysis_mode:

  Single or pairwise analysis mode.

- contrast:

  Contrast vector for single analysis.

- design_formula:

  Design formula.

- pairwise_group_cols:

  Metadata columns used to build pairwise groups.

- pairwise_contrast_col:

  Name of pairwise contrast column.

- pairwise_contrasts:

  Optional named list of pairwise contrasts.

- filter_method:

  Gene filtering method.

- pairwise_mode:

  Pairwise contrast generation mode.

- min_count:

  Minimum count threshold.

- min_samples:

  Minimum number of samples passing `min_count`.

- min_total:

  Minimum total count.

- rhythmicity_analysis:

  Logical. Run MetaCycle + cosinor rhythmicity analysis on the fitted
  `dds` object (see
  [`run_deggo_rhythmicity()`](https://ymbouamboua.github.io/DEGgo/reference/run_deggo_rhythmicity.md)).
  Disabled by default; requires a numeric time column in `metadata`.

- rhythmicity_time_col:

  Metadata column with numeric time (e.g. ZT hours) used for rhythmicity
  analysis.

- rhythmicity_group_col:

  Optional metadata column defining exactly two groups for the
  differential rhythmicity test.

- rhythmicity_assay:

  Assay extracted from `dds` for rhythmicity analysis. One of `"vst"`,
  `"normalized"`, `"log2_normalized"`, `"raw"`.

- rhythmicity_methods:

  Rhythmicity methods to run: `"meta2d"`, `"cosinor"`, or both.

- rhythmicity_period_range:

  MetaCycle period search window.

- rhythmicity_cycle_length:

  Assumed period for the cosinor fit.

- rhythmicity_cycMethod:

  MetaCycle methods to combine.

- cosinor_engine:

  Cosinor fitting engine: `"auto"`, `"package"` (`cosinor`/`cosinor2`),
  or `"manual"` (dependency-free
  [`lm()`](https://rdrr.io/r/stats/lm.html) fallback).

- rhythmicity_plots:

  Logical. Generate rhythmicity diagnostic plots.

- rhythmicity_n_top_plots:

  Number of top rhythmic genes to plot.

- generate_report:

  Logical. Generate report.

- report_formats:

  Report formats.

- report_template:

  Optional report template path.

- generate_pptx:

  Logical. Generate PowerPoint.

- pptx_file:

  Optional PowerPoint output file.

- save_reproducibility:

  Logical. Save reproducibility bundle.

- save_clean_inputs:

  Logical. Save cleaned input tables.

- txtsize:

  Base text size.

- heatmap_annotation_cols:

  Character vector of metadata columns used as heatmap annotations, or
  `"auto"` for automatic selection.

- palette:

  Optional color palette used by DEGgo plots. Can be a named DEGgo
  palette or a character vector of colors.

- seed:

  Random seed.

## Value

A DEGgo results object.
