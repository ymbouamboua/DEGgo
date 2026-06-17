DEGgo
================

- [DEGgo](#deggo)
  - [Overview](#overview)
  - [Installation](#installation)
  - [Workflow](#workflow)
  - [Input data](#input-data)
    - [Count table](#count-table)
    - [Metadata](#metadata)
  - [Main functions](#main-functions)
- [Complete DEGgo workflow](#complete-deggo-workflow)
  - [1. Raw count quality control](#1-raw-count-quality-control)
  - [2. Remove flagged samples](#2-remove-flagged-samples)
  - [3. Save cleaned tables](#3-save-cleaned-tables)
  - [4. Re-run QC after cleaning](#4-re-run-qc-after-cleaning)
  - [5. Tissue marker validation](#5-tissue-marker-validation)
- [Differential expression analysis](#differential-expression-analysis)
  - [6A. Single comparison mode](#6a-single-comparison-mode)
  - [6B. Pairwise contrast mode](#6b-pairwise-contrast-mode)
- [Automated HTML/PDF report](#automated-htmlpdf-report)
  - [HTML report](#html-report)
  - [HTML and PDF report](#html-and-pdf-report)
  - [Generate a report from existing
    results](#generate-a-report-from-existing-results)
- [Gene expression extraction and
  plotting](#gene-expression-extraction-and-plotting)
  - [Extract expression values](#extract-expression-values)
  - [Plot one gene](#plot-one-gene)
  - [Plot another gene](#plot-another-gene)
- [GO enrichment](#go-enrichment)
  - [Plot GO terms](#plot-go-terms)
  - [Plot all GO terms across pairwise
    comparisons](#plot-all-go-terms-across-pairwise-comparisons)
- [Extract genes and GO terms of
  interest](#extract-genes-and-go-terms-of-interest)
  - [1. Define genes of interest](#1-define-genes-of-interest)
  - [2. Extract DEG genes](#2-extract-deg-genes)
  - [3. Extract GO terms containing genes of
    interest](#3-extract-go-terms-containing-genes-of-interest)
  - [4. Extract GO terms using biological
    keywords](#4-extract-go-terms-using-biological-keywords)
- [Output structure](#output-structure)
- [Citation](#citation)
- [License](#license)

<p align="center">

<img src="man/figures/DEGgo_logo.png" width="300">
</p>

# DEGgo

<p align="center">

<img src="https://img.shields.io/badge/R-%3E%3D4.3-blue">
<img src="https://img.shields.io/badge/license-MIT-green">
<img src="https://img.shields.io/github/stars/ymbouamboua/DEGgo?style=social">
</p>

## Overview

**DEGgo** is an R package for end-to-end bulk RNA-seq differential
expression analysis and Gene Ontology enrichment.

DEGgo provides a complete workflow from raw count matrices to biological
interpretation:

- raw count quality control;
- sample validation and filtering;
- tissue marker verification;
- differential expression analysis;
- single and pairwise contrasts;
- gene ontology enrichment;
- extraction of genes and GO terms of interest;
- publication-ready visualization;
- automated HTML/PDF reporting.

## Installation

``` r
install.packages("remotes")
remotes::install_github("ymbouamboua/DEGgo")
```

``` r
library(DEGgo)
```

## Workflow

``` text
Raw counts + metadata
        │
        ▼
Input preparation and sample matching
        │
        ▼
Quality control and marker validation
        │
        ▼
Low-expression filtering
        │
        ▼
Differential expression analysis
        │
        ├── Single contrast mode
        │
        └── Pairwise contrast mode
                 │
                 ▼
Annotation and DEG export
                 │
                 ▼
Volcano plots, PCA and heatmaps
                 │
                 ▼
Gene Ontology enrichment
                 │
                 ▼
Gene/GO extraction and visualization
                 │
                 ▼
HTML/PDF DEGgo report
```

## Input data

### Count table

DEGgo accepts raw count tables or matrices. Gene annotation columns can
be automatically detected using `gene_col` and `feature_col`.

``` text
gene_id      gene_name    Sample1    Sample2    Sample3
ENSMUSG1     Adipoq       120        145        98
ENSMUSG2     Lep          65         80         50
ENSMUSG3     Ucp1         12         18         250
```

### Metadata

Metadata must contain one row per sample. A `sample` column is expected
or automatically detected using `sample_col`.

``` text
sample      tissue    sex       treatment
Sample1     WAT       Female    PBS
Sample2     WAT       Female    PAMH
Sample3     BAT       Female    PBS
```

For pairwise analysis, DEGgo can build a `condition` variable from
multiple metadata columns, for example:

``` r
pairwise_group_cols = c("treatment", "sex", "tissue")
```

This creates groups such as:

``` text
PBS_Female_WAT
PAMH_Female_WAT
PBS_Male_BAT
PAMH_Male_BAT
```

## Main functions

| Function | Description |
|----|----|
| `explore_bulk_rnaseq()` | Raw and cleaned RNA-seq quality control |
| `remove_flagged_samples()` | Remove failed or low-quality samples |
| `marker_score_check()` | Tissue marker scoring and sample swap detection |
| `plot_gene_heatmap()` | Marker or selected-gene heatmap |
| `run_deggo()` | Main DEGgo differential expression workflow |
| `run_go_enrichment()` | Gene Ontology enrichment for DEG results |
| `plot_go_terms()` | Publication-ready GO plot |
| `plot_all_go_terms()` | Plot GO terms across multiple comparisons |
| `extract_expression()` | Extract raw, normalized, log2-normalized or VST expression values |
| `plot_gene_expression()` | Plot normalized gene expression for selected genes |
| `generate_deggo_report()` | Generate HTML/PDF DEGgo report |
| `deggo_extract_deg_genes()` | Extract genes of interest from DEG results |
| `deggo_extract_go_genes_pairwise()` | Extract GO terms containing selected genes |
| `deggo_extract_go_keywords()` | Extract GO terms matching biological keywords |

# Complete DEGgo workflow

## 1. Raw count quality control

``` r
qc <- explore_bulk_rnaseq(
  counts = counts,
  metadata = metadata,
  markers = c(
    "Adipoq", "Lep", "Ucp1",
    "Amh", "Amhr2", "Ar", "Akr1c3"
  ),
  output_dir = file.path(outdir, "DEGgo_QC_raw")
)

qc$qc
```

## 2. Remove flagged samples

``` r
cleaned <- remove_flagged_samples(
  counts = counts,
  metadata = metadata,
  qc_table = qc$qc,
  sample_col = "sample",
  remove_col = "recommend_remove",
  gene_cols = c("gene_id", "gene_name"),
  verbose = TRUE
)

counts <- cleaned$counts
metadata <- cleaned$metadata
```

## 3. Save cleaned tables

``` r
savetbl(
  counts,
  file.path(root, "data", "processed", "clean_counts"),
  "xlsx"
)

savetbl(
  metadata,
  file.path(root, "data", "processed", "clean_metadata"),
  "xlsx"
)
```

## 4. Re-run QC after cleaning

``` r
qc_clean <- explore_bulk_rnaseq(
  counts = counts,
  metadata = metadata,
  markers = c(
    "Adipoq", "Lep", "Ucp1",
    "Amh", "Amhr2", "Ar", "Akr1c3"
  ),
  output_dir = file.path(outdir, "DEGgo_QC_clean")
)

qc_clean$qc
```

## 5. Tissue marker validation

``` r
marker_sets <- list(
  BAT = c(
    "Ucp1", "Cidea", "Ppargc1a",
    "Prdm16", "Dio2", "Elovl3"
  ),
  WAT = c(
    "Adipoq", "Lep", "Retn",
    "Fabp4", "Pparg", "Lpl"
  ),
  TESTIS = c(
    "Amh", "Amhr2", "Sox9", "Dhh",
    "Star", "Cyp11a1", "Cyp17a1",
    "Insl3", "Gata4", "Oct4", "Ddx4"
  ),
  OVARY = c(
    "Foxl2", "Fshr", "Cyp19a1",
    "Bmp15", "Gdf9"
  )
)

marker_genes <- unique(unlist(marker_sets))

plot_gene_heatmap(
  counts = counts,
  metadata = metadata,
  genes = marker_genes,
  gene_col = "gene_id",
  feature_col = "gene_name",
  sample_col = "sample",
  annotation_cols = c("tissue", "sex", "treatment"),
  order_by = c("tissue", "sex", "treatment"),
  output_dir = file.path(outdir, "DEGgo_QC_clean"),
  filename = "Tissue_Marker_Check_Heatmap",
  main = "Tissue marker check"
)

marker_check <- marker_score_check(
  counts = counts,
  metadata = metadata,
  marker_sets = marker_sets,
  sample_col = "sample",
  group_col = "tissue",
  feature_col = "gene_name"
)

marker_check$scores
marker_check$swaps
marker_check$plot
```

# Differential expression analysis

## 6A. Single comparison mode

Use `analysis_mode = "single"` when you have one contrast such as `PAMH`
vs `PBS`, `KO` vs `WT`, or `treated` vs `control`.

The metadata must contain a `condition` column when using single mode.

``` r
results_single <- run_deggo(
  counts = counts,
  metadata = metadata,
  organism = "mouse",
  method = "DESeq2",
  analysis_mode = "single",
  contrast = c("condition", "PAMH", "PBS"),
  design_formula = ~ condition,
  filter_method = "count",
  min_count = 2,
  min_samples = 2,
  min_total = 10,
  generate_report = TRUE,
  report_formats = "html",
  output_dir = file.path(outdir, "DEGgo_single")
)

results_single$summary
results_single$report_files
```

## 6B. Pairwise contrast mode

Use `analysis_mode = "pairwise"` when your experiment contains multiple
structured comparisons, such as treatment by sex by tissue.

``` r
pairwise_contrasts <- list(
  WAT_Female_PAMH_vs_PBS =
    c("comparison_group", "PAMH_Female_WAT", "PBS_Female_WAT"),

  WAT_Male_PAMH_vs_PBS =
    c("comparison_group", "PAMH_Male_WAT", "PBS_Male_WAT"),

  BAT_Female_PAMH_vs_PBS =
    c("comparison_group", "PAMH_Female_BAT", "PBS_Female_BAT"),

  BAT_Male_PAMH_vs_PBS =
    c("comparison_group", "PAMH_Male_BAT", "PBS_Male_BAT"),

  Ovary_PAMH_vs_PBS =
    c("comparison_group", "PAMH_Female_OVARY", "PBS_Female_OVARY"),

  Testis_PAMH_vs_PBS =
    c("comparison_group", "PAMH_Male_TESTIS", "PBS_Male_TESTIS")
)
```

``` r
results <- run_deggo(
  counts = counts,
  metadata = metadata,
  organism = "mouse",
  method = "DESeq2",
  analysis_mode = "pairwise",
  pairwise_group_cols = c("treatment", "sex", "tissue"),
  pairwise_contrast_col = "comparison_group",
  pairwise_contrasts = pairwise_contrasts,
  filter_method = "count",
  min_count = 2,
  min_samples = 2,
  min_total = 10,
  generate_report = TRUE,
  report_formats = "html",
  output_dir = outdir
)

results$summary
results$report_files
```

``` r
saveRDS(
  results,
  file.path(outdir, "DEGgo_results.rds")
)
```

# Automated HTML/PDF report

DEGgo can automatically generate an HTML or PDF report at the end of
`run_deggo()` using the arguments `generate_report`, `report_formats`,
and `report_template`.

## HTML report

``` r
results <- run_deggo(
  counts = counts,
  metadata = metadata,
  organism = "mouse",
  method = "DESeq2",
  analysis_mode = "pairwise",
  pairwise_group_cols = c("treatment", "sex", "tissue"),
  pairwise_contrast_col = "comparison_group",
  pairwise_contrasts = pairwise_contrasts,
  generate_report = TRUE,
  report_formats = "html",
  output_dir = outdir
)

results$report_files$html
```

## HTML and PDF report

``` r
results <- run_deggo(
  counts = counts,
  metadata = metadata,
  organism = "mouse",
  method = "DESeq2",
  analysis_mode = "pairwise",
  pairwise_group_cols = c("treatment", "sex", "tissue"),
  pairwise_contrast_col = "comparison_group",
  pairwise_contrasts = pairwise_contrasts,
  generate_report = TRUE,
  report_formats = c("html", "pdf"),
  output_dir = outdir
)

results$report_files
```

PDF generation requires a LaTeX installation, for example TinyTeX.

``` r
install.packages("tinytex")
tinytex::install_tinytex()
```

## Generate a report from existing results

If you already have a DEGgo results object, you can generate the report
directly:

``` r
report_files <- generate_deggo_report(
  results = results,
  output_dir = results$output_dir,
  formats = c("html", "pdf")
)

report_files
```

# Gene expression extraction and plotting

## Extract expression values

`extract_expression()` extracts raw, normalized, log2-normalized, or
variance-stabilized expression values from a DESeq2 object.

``` r
expr_df <- extract_expression(
  dds = results$dds,
  metadata = results$metadata,
  genes = c("Adipoq", "Lep", "Ucp1", "Amh", "Amhr2"),
  assay = "vst",
  gene_col = "SYMBOL"
)

head(expr_df)
```

## Plot one gene

``` r
plot_gene_expression(
  expr_df,
  gene = "Adipoq",
  x = "treatment",
  color = "treatment",
  facet = "tissue",
  geom = "violin"
)
```

## Plot another gene

``` r
plot_gene_expression(
  expr_df,
  gene = "Ucp1",
  x = "treatment",
  color = "treatment",
  facet = "tissue",
  geom = "boxplot"
)
```

# GO enrichment

DEGgo performs Gene Ontology enrichment on significant differentially
expressed genes.

``` r
go <- run_go_enrichment(
  sig_deg = results$sig_deg$WAT_Female_PAMH_vs_PBS,
  comparison = "WAT_Female_PAMH_vs_PBS",
  ontology = "BP",
  orgdb = org.Mm.eg.db::org.Mm.eg.db,
  output_dir = file.path(outdir, "GO_WAT_Female_PAMH_vs_PBS")
)

go$go_results
go$go_plot
```

## Plot GO terms

``` r
plot_go_terms(
  go_df = go$go_results,
  comparison = "WAT Female PAMH vs PBS",
  top_n = 15
)
```

## Plot all GO terms across pairwise comparisons

``` r
plot_all_go_terms(
  results = results,
  top_n = 10,
  output_dir = file.path(outdir, "pairwise_GO_plots")
)
```

# Extract genes and GO terms of interest

## 1. Define genes of interest

``` r
genes_interest <- c(
  "Amh", "Amhr2",
  "Ar", "Akr1c3",
  "Star", "Cyp11a1", "Cyp17a1", "Cyp19a1",
  "Foxl2", "Fshr",
  "Sox9", "Insl3",
  "Pparg", "Fabp4", "Adipoq", "Lep",
  "Ucp1", "Prdm16", "Ppargc1a", "Cidea"
)
```

## 2. Extract DEG genes

``` r
gene_summary <- deggo_extract_deg_genes(
  results,
  genes_interest
)

gene_summary <- gene_summary[
  order(gene_summary$comparison, gene_summary$SYMBOL),
]

gene_summary_sig <- subset(
  gene_summary,
  significance %in% c("Up", "Down")
)

gene_summary
gene_summary_sig
```

## 3. Extract GO terms containing genes of interest

``` r
go_genes_pairwise <- deggo_extract_go_genes_pairwise(
  results,
  genes_interest
)

go_genes_pairwise
```

## 4. Extract GO terms using biological keywords

``` r
go_keywords <- c(
  "tgf", "bmp", "smad",
  "inflammation", "immune",
  "steroidogenesis", "androgen", "steroid",
  "epigenetic",
  "adipocyte", "lipid", "fatty acid",
  "thermogenesis", "mitochond",
  "hormone", "reproductive",
  "gonad", "ovary", "testis"
)

go_keywords_results <- deggo_extract_go_keywords(
  results,
  keywords = go_keywords
)

go_keywords_results
```

# Output structure

For pairwise DESeq2 analysis, DEGgo creates a dated analysis directory:

``` text
DEGgo_results/
└── DEGgo_YYYY_MM_DD_DESeq2_pairwise/
    │
    ├── pairwise_DESeq2_results/
    ├── pairwise_DESeq2_significant/
    ├── pairwise_volcano/
    ├── pairwise_heatmaps/
    ├── pairwise_PCA/
    ├── pairwise_GO/
    ├── pairwise_GO_plots/
    ├── DEGgo_Report.html
    ├── DEGgo_Report.pdf
    ├── pairwise_summary.tsv
    ├── session_info.txt
    └── DEGgo_results.rds
```

For single analysis, the output directory has the same general structure
but contains single-comparison result files.

# Citation

If you use DEGgo in your work, please cite:

> Yvon MBOUAMBOUA and Paolo Giacobini. DEGgo: an integrated framework
> for bulk RNA-seq differential expression analysis and functional
> enrichment.

# License

MIT © Yvon MBOUAMBOUA
