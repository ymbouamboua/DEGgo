DEGgo
================

<p align="center">

<img src="man/figures/DEGgo_logo.svg" alt="DEGgo logo" width="500">
</p>

<p align="center">

An integrated framework for automated bulk RNA-seq differential
expression analysis, functional enrichment, visualization, and
reproducible reporting.
</p>

<p align="center">

<img src="https://img.shields.io/badge/R-%3E%3D4.3-blue">
<img src="https://img.shields.io/badge/license-MIT-green"> </a>
<a href="https://doi.org/10.5281/zenodo.20785178">
<img src="https://zenodo.org/badge/DOI/10.5281/zenodo.20785178.svg">
<img src="https://img.shields.io/github/stars/ymbouamboua/DEGgo?style=social">
</a>
</p>

## Overview

**DEGgo** is an R package for automated bulk RNA-seq downstream
analysis.

Starting from a raw count matrix and sample metadata, DEGgo performs
sample quality control, differential expression analysis using
**DESeq2**, **edgeR**, or **limma**, Gene Ontology enrichment, circadian
rhythmicity analysis, publication-ready visualization, and automated
HTML, PDF, and PowerPoint reporting within a reproducible workflow.

For time-course experiments, `run_deggo_rhythmicity()` provides
rhythmicity detection using **MetaCycle** and single-component **cosinor
regression**. The function can also test differential rhythmicity
between two biological groups when the experimental design supports such
a comparison.

The package is designed for bioinformaticians and experimental
biologists, minimizing repetitive scripting while preserving analytical
flexibility.

## Features

- Automated bulk RNA-seq workflow
- Raw and post-filtering sample quality control
- Differential expression with DESeq2, edgeR, or limma
- Single-comparison and pairwise analyses
- Gene Ontology over-representation analysis
- Circadian and time-course rhythmicity analysis
- MetaCycle rhythmicity detection
- Single-component cosinor regression
- Optional differential-rhythmicity testing between two groups
- Gene-symbol annotation of rhythmicity results
- Publication-ready PCA, volcano plots, heatmaps, and rhythmicity plots
- Automated HTML, PDF, and PowerPoint reports
- Reproducibility files, software versions, and session information
- Human, mouse, rat, and custom `OrgDb` support

## Workflow

<p align="center">

<img src="man/figures/DEGgo_workflow.png" width="850">
</p>

The main DEGgo workflow includes:

1.  Input validation and matching of counts and metadata
2.  Raw sample quality control
3.  Count filtering and optional removal of flagged samples
4.  Differential expression analysis
5.  PCA, heatmap, and volcano plot generation
6.  Gene Ontology enrichment
7.  Automated report and PowerPoint generation
8.  Reproducibility and session-information export

For time-course experiments, the independent rhythmicity module adds:

1.  Time-point validation and sample matching
2.  Expression-matrix preparation
3.  MetaCycle rhythmicity detection
4.  Cosinor model fitting
5.  Optional differential-rhythmicity testing
6.  Consensus rhythmic-gene classification
7.  Rhythmicity tables and diagnostic plots

## Installation

``` r
install.packages("remotes")

# Development version
remotes::install_github("ymbouamboua/DEGgo")
```

``` r
library(DEGgo)
```

Optional rhythmicity dependencies can be installed with:

``` r
install.packages("MetaCycle")

# Optional cosinor engines
install.packages("cosinor")
install.packages("cosinor2")
```

When the optional cosinor packages are unavailable,
`run_deggo_rhythmicity()` can use its internal `lm()`-based
implementation.

## Quick start

The package includes a ready-to-use RNA-seq dataset derived from the
Bioconductor **airway** experiment.

``` r
counts <- read.delim(
  system.file(
    "extdata",
    "airway_counts.tsv",
    package = "DEGgo"
  ),
  check.names = FALSE
)

metadata <- read.delim(
  system.file(
    "extdata",
    "airway_metadata.tsv",
    package = "DEGgo"
  ),
  check.names = FALSE
)
```

## Single differential expression analysis

``` r
results <- run_deggo(
  counts = counts,
  metadata = metadata,
  organism = "human",
  method = "DESeq2",
  analysis_mode = "single",
  sample_col = "sample",
  design_formula = ~ treatment,
  contrast = c(
    "treatment",
    "treated",
    "untreated"
  ),
  output_dir = "DEGgo_results"
)

results$summary
```

The comparison is interpreted as:

``` text
treated versus untreated
```

Therefore, positive `log2FoldChange` values represent higher expression
in the treated group.

## Pairwise differential expression analysis

Create a combined comparison variable when pairwise comparisons involve
several metadata columns:

``` r
metadata$group <- interaction(
  metadata$condition,
  metadata$sex,
  sep = "_",
  drop = TRUE
)
```

Define the comparisons:

``` r
pairwise_contrasts <- list(
  TreatmentA_vs_Control_Male = c(
    "group",
    "TreatmentA_Male",
    "Control_Male"
  ),

  TreatmentA_vs_Control_Female = c(
    "group",
    "TreatmentA_Female",
    "Control_Female"
  ),

  Male_vs_Female_Control = c(
    "group",
    "Control_Male",
    "Control_Female"
  )
)
```

Run the pairwise workflow:

``` r
pairwise_results <- run_deggo(
  counts = counts,
  metadata = metadata,
  organism = "mouse",
  method = "DESeq2",
  analysis_mode = "pairwise",
  sample_col = "sample",
  design_formula = ~ group,
  pairwise_group_cols = c(
    "condition",
    "sex"
  ),
  pairwise_contrast_col = "group",
  pairwise_contrasts = pairwise_contrasts,
  output_dir = "DEGgo_pairwise_results"
)
```

## Circadian rhythmicity analysis

`run_deggo_rhythmicity()` is a standalone workflow for detecting
rhythmic genes in circadian or periodic time-course experiments.

It accepts:

- a numeric genes-by-samples matrix;
- a data frame containing expression values; or
- a DESeq2 `DESeqDataSet`.

The metadata must contain:

- a sample identifier column;
- a numeric time column;
- optionally, a two-level group column for differential-rhythmicity
  testing.

### Prepare circadian metadata

Suppose the metadata contain Zeitgeber Time values such as `ZT2`, `ZT6`,
`ZT10`, `ZT14`, `ZT18`, and `ZT22`.

``` r
rhythm_metadata <- metadata

rhythm_metadata$time <- as.numeric(
  sub(
    "^ZT",
    "",
    rhythm_metadata$ZT
  )
)

rhythm_metadata$time
```

Expected values:

``` text
2 6 10 14 18 22
```

The time variable must remain numeric.

### Run rhythmicity analysis from a DESeq2 object

A DESeq2 object returned by `run_deggo()` can be reused directly:

``` r
rhythm_results <- run_deggo_rhythmicity(
  expr = pairwise_results$dds,
  metadata = rhythm_metadata,
  sample_col = "sample",
  time_col = "time",
  group_col = NULL,
  assay = "vst",
  methods = c(
    "meta2d",
    "cosinor"
  ),
  period_range = c(
    20,
    28
  ),
  cycle_length = 24,
  cycMethod = c(
    "ARS",
    "JTK",
    "LS"
  ),
  padj_cutoff = 0.05,
  cosinor_engine = "auto",
  output_dir = "DEGgo_rhythmicity_results",
  project_name = "Circadian RNA-seq analysis",
  generate_plots = TRUE,
  n_top_plots = 25,
  seed = 4173
)
```

### Add gene-symbol annotation

If the expression matrix uses Ensembl identifiers, provide a
gene-annotation table so that result tables and plots display readable
symbols:

``` r
gene_annotation <- unique(
  counts[
    ,
    c(
      "gene_id",
      "Gene"
    )
  ]
)

rhythm_results <- run_deggo_rhythmicity(
  expr = pairwise_results$dds,
  metadata = rhythm_metadata,
  sample_col = "sample",
  time_col = "time",
  assay = "vst",
  methods = c(
    "meta2d",
    "cosinor"
  ),
  period_range = c(
    20,
    28
  ),
  cycle_length = 24,
  padj_cutoff = 0.05,
  cosinor_engine = "auto",
  gene_annotation = gene_annotation,
  gene_id_col = "gene_id",
  gene_symbol_col = "Gene",
  show_gene_id = FALSE,
  output_dir = "DEGgo_rhythmicity_results",
  generate_plots = TRUE,
  n_top_plots = 25
)
```

The annotation table is supplied by the analysis script. The package
function does not depend on an external global object named `counts`.

### Inspect rhythmicity results

The combined result table is available through:

``` r
head(rhythm_results$summary)
```

Method-specific tables are available through:

``` r
head(rhythm_results$metacycle)
head(rhythm_results$cosinor)
```

Summary statistics are stored in:

``` r
rhythm_results$stats
```

The result object contains:

``` r
names(rhythm_results)
```

Typical components include:

``` text
summary
metacycle
cosinor
cosinor_diff_group
metadata
expr
gene_annotation
stats
params
output_dir
output_dirs
plots
version
```

### Extract consensus rhythmic genes

Genes detected as rhythmic by both MetaCycle and cosinor can be
extracted from the combined summary:

``` r
consensus_rhythmic <- subset(
  rhythm_results$summary,
  rhythmic_by == "both"
)

consensus_rhythmic <- consensus_rhythmic[
  order(
    consensus_rhythmic$cosinor_padj,
    consensus_rhythmic$metacycle_padj
  ),
]

head(consensus_rhythmic)
```

Save the consensus table:

``` r
write.table(
  consensus_rhythmic,
  file = file.path(
    rhythm_results$output_dirs$results,
    "consensus_rhythmic_genes.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
```

### Differential rhythmicity between two groups

Differential rhythmicity can be tested when `group_col` contains exactly
two biologically independent groups.

For example:

``` r
rhythm_results_grouped <- run_deggo_rhythmicity(
  expr = pairwise_results$dds,
  metadata = rhythm_metadata,
  sample_col = "sample",
  time_col = "time",
  group_col = "condition",
  assay = "vst",
  methods = c(
    "meta2d",
    "cosinor"
  ),
  period_range = c(
    20,
    28
  ),
  cycle_length = 24,
  padj_cutoff = 0.05,
  cosinor_engine = "auto",
  output_dir = "DEGgo_differential_rhythmicity"
)
```

The differential-rhythmicity table is available through:

``` r
head(
  rhythm_results_grouped$cosinor_diff_group
)
```

A group variable should not be derived directly from the same time
variable used in the cosinor model. For example, dividing circadian time
points into `Light` and `Dark` categories usually creates strong
confounding between group and time and is not an appropriate
differential-rhythmicity design.

A valid design would instead compare two biological conditions, each
sampled across the complete circadian cycle.

### Rhythmicity output files

The rhythmicity workflow produces files such as:

``` text
DEGgo_rhythmicity_results/
├── plots/
│   ├── rhythmicity_period_distribution.png
│   ├── rhythmicity_pvalue_comparison.png
│   └── individual_genes/
├── results/
│   ├── metacycle_results.tsv
│   ├── cosinor_results.tsv
│   ├── cosinor_differential_rhythmicity.tsv
│   └── rhythmicity_summary.tsv
└── reproducibility/
```

The exact files depend on the selected methods and whether grouped
differential-rhythmicity analysis was requested.

## Main functions

| Function | Description |
|:---|:---|
| `run_deggo()` | Run the complete bulk RNA-seq differential expression workflow |
| `run_deggo_rhythmicity()` | Detect rhythmic genes using MetaCycle and cosinor regression |
| `check_raw_counts()` | Validate raw count matrices and sample identifiers before analysis |
| `explore_bulk_rnaseq()` | Perform exploratory quality control on raw RNA-seq data |
| `remove_flagged_samples()` | Remove low-quality samples identified during quality control |
| `marker_score_check()` | Assess marker gene expression for biological sample validation |
| `run_sample_qc()` | Perform post-filtering sample quality control and diagnostics |
| `run_go_enrichment()` | Perform Gene Ontology over-representation analysis |
| `plot_go_terms()` | Generate publication-ready Gene Ontology enrichment plots |
| `plot_all_go_terms()` | Compare Gene Ontology enrichment across multiple analyses |
| `extract_expression()` | Extract raw, normalized, log2-normalized, or VST expression values |
| `plot_gene_expression()` | Visualize selected genes across biological groups |
| `plot_gene_heatmap()` | Visualize marker or selected genes as a heatmap |
| `plot_heatmap()` | Generate heatmaps of differentially expressed genes |
| `plot_pca()` | Perform principal component analysis |
| `plot_volcano()` | Generate publication-ready volcano plots |
| `generate_deggo_report()` | Generate an automated HTML or PDF analysis report |
| `generate_deggo_pptx()` | Generate a PowerPoint summary report |
| `deggo_extract_deg_genes()` | Extract selected genes from differential expression results |
| `deggo_extract_go_genes()` | Extract genes associated with enriched Gene Ontology terms |
| `deggo_extract_go_genes_pairwise()` | Extract GO-associated genes across pairwise comparisons |
| `deggo_extract_go_keywords()` | Search enriched Gene Ontology terms using biological keywords |

## Rhythmicity function

The principal rhythmicity function has the following interface:

``` r
run_deggo_rhythmicity(
  expr,
  metadata,
  sample_col = "sample",
  time_col = "time",
  group_col = NULL,
  assay = c(
    "vst",
    "normalized",
    "log2_normalized",
    "raw"
  ),
  methods = c(
    "meta2d",
    "cosinor"
  ),
  period_range = c(
    20,
    28
  ),
  cycle_length = 24,
  cycMethod = c(
    "ARS",
    "JTK",
    "LS"
  ),
  padj_cutoff = 0.05,
  cosinor_engine = c(
    "auto",
    "package",
    "manual"
  ),
  gene_annotation = NULL,
  gene_id_col = "gene_id",
  gene_symbol_col = NULL,
  show_gene_id = FALSE,
  output_dir = "DEGgo_rhythmicity_out",
  project_name = NULL,
  generate_plots = TRUE,
  n_top_plots = 20,
  seed = 4173,
  verbose = TRUE
)
```

### Important rhythmicity parameters

| Parameter | Description |
|:---|:---|
| `expr` | Numeric genes-by-samples matrix, data frame, or DESeq2 object |
| `metadata` | Sample metadata containing sample identifiers and numeric time |
| `sample_col` | Metadata column containing sample names |
| `time_col` | Numeric time variable |
| `group_col` | Optional two-level biological grouping variable |
| `assay` | Expression transformation used for a DESeq2 object |
| `methods` | Rhythmicity methods: `meta2d`, `cosinor`, or both |
| `period_range` | Minimum and maximum periods tested by MetaCycle |
| `cycle_length` | Fixed period used in the cosinor model |
| `cycMethod` | MetaCycle algorithms to combine |
| `padj_cutoff` | Adjusted p-value threshold |
| `cosinor_engine` | Package-based, manual, or automatically selected engine |
| `gene_annotation` | Optional gene identifier-to-symbol annotation table |
| `gene_id_col` | Gene identifier column in the annotation table |
| `gene_symbol_col` | Gene-symbol column in the annotation table |
| `show_gene_id` | Include original identifiers in rhythmicity plot titles |
| `generate_plots` | Generate diagnostic and individual-gene plots |
| `n_top_plots` | Number of top rhythmic genes to plot |

## Supported organisms

DEGgo provides built-in annotation support for the following organisms:

| Organism                  | Parameter | OrgDb package  |
|:--------------------------|:----------|:---------------|
| Human (*Homo sapiens*)    | `"human"` | `org.Hs.eg.db` |
| Mouse (*Mus musculus*)    | `"mouse"` | `org.Mm.eg.db` |
| Rat (*Rattus norvegicus*) | `"rat"`   | `org.Rn.eg.db` |

Additional organisms can be analyzed by supplying a compatible
Bioconductor `OrgDb` annotation database.

Rhythmicity analysis itself does not require an `OrgDb` package when the
expression matrix already contains suitable gene identifiers. A custom
annotation table can be provided through `gene_annotation`.

## Documentation

A complete tutorial and advanced examples are available in the package
vignette.

``` r
browseVignettes("DEGgo")
```

Function-specific documentation can be opened with:

``` r
?run_deggo
?run_deggo_rhythmicity
?plot_volcano
?plot_heatmap
```

## Citation

If you use **DEGgo** in your research, please cite:

> Yvon Mbouamboua, Vincent Prevot, Paolo Giacobini (2026). **DEGgo:
> automated bulk RNA-seq differential expression, functional enrichment,
> and rhythmicity analysis.** Zenodo.
> <https://doi.org/10.5281/zenodo.20785178>

## License

MIT © Yvon MBOUAMBOUA
