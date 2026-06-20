# Run DEGgo bulk RNA-seq downstream analysis

Run an automated bulk RNA-seq downstream analysis workflow including
input preparation, sample matching, differential expression analysis,
result annotation, visualization, heatmap generation, PCA analysis, Gene
Ontology enrichment, and optional HTML/PDF report generation.

## Usage

``` r
run_deggo(
  counts,
  metadata,
  gene_col = c("gene_id", "GeneID", "gene", "Gene", "ENSEMBL", "ensembl", "ensembl_id"),
  feature_col = c("gene_name", "SYMBOL", "symbol", "gene_symbol", "external_gene_name"),
  sample_col = c("sample", "Sample", "SAMPLE"),
  prepare_input = TRUE,
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
  pairwise_auto = TRUE,
  pairwise_mode = c("all", "within_first", "within_second"),
  min_count = 5,
  min_samples = 2,
  min_total = 10,
  generate_report = TRUE,
  report_formats = "html",
  report_template = NULL,
  save_reproducibility = TRUE,
  save_clean_inputs = TRUE,
  seed = 123
)
```

## Arguments

- counts:

  Raw count table or matrix.

- metadata:

  Sample metadata data frame.

- gene_col:

  Character vector of possible gene identifier columns.

- feature_col:

  Character vector of possible feature annotation columns.

- sample_col:

  Character vector of possible sample identifier columns.

- prepare_input:

  Logical. If TRUE, automatically matches counts and metadata.

- output_dir:

  Output directory. If NULL, a dated directory is created.

- padj_cutoff:

  Adjusted p-value cutoff.

- logfc_cutoff:

  Absolute log2 fold-change cutoff.

- top_n_heatmap:

  Number of top genes for heatmap.

- top_n_labels:

  Number of top genes labeled in volcano plot.

- ontology:

  GO ontology: BP, MF, or CC.

- organism:

  Organism: human, mouse, rat or custom by user.

- orgdb:

  Optional OrgDb object.

- method:

  DE method: DESeq2, edgeR, or limma.

- analysis_mode:

  Analysis mode: single or pairwise.

- contrast:

  Optional contrast for single analysis.

- design_formula:

  Design formula.

- pairwise_group_cols:

  Metadata columns for pairwise groups.

- pairwise_contrast_col:

  Name of pairwise contrast column.

- pairwise_contrasts:

  Optional named list of contrasts.

- filter_method:

  Filtering method: count, cpm, or none.

- pairwise_auto:

  Automatically generate pairwise contrasts.

- pairwise_mode:

  Pairwise mode: all, within_first, or within_second.

- min_count:

  Minimum count for filtering.

- min_samples:

  Minimum samples for filtering.

- min_total:

  Minimum total count for filtering.

- generate_report:

  Logical. Generate DEGgo report.

- report_formats:

  Character vector. Report formats: html and/or pdf.

- report_template:

  Optional path to report R Markdown template.

- save_reproducibility:

  Logical. Save run parameters and reproducibility files.

- save_clean_inputs:

  Logical. Save cleaned counts and matched metadata.

- seed:

  Random seed.

## Value

A list containing DEG results, plots, GO results, reports, and output
paths.
