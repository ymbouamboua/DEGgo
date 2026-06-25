# DEGgo

![DEGgo logo](reference/figures/DEGgo_logo.png)

An integrated framework for automated bulk RNA-seq differential
expression analysis, functional enrichment, visualization, and
reproducible reporting.

[![](https://zenodo.org/badge/DOI/10.5281/zenodo.20785178.svg)](https://doi.org/10.5281/zenodo.20785178)
![](https://img.shields.io/badge/R-%3E%3D4.3-blue)![](https://img.shields.io/badge/license-MIT-green)![](https://img.shields.io/github/stars/ymbouamboua/DEGgo?style=social)

## Overview

**DEGgo** is an R package for automated bulk RNA-seq downstream
analysis.

Starting from a raw count matrix and sample metadata, DEGgo performs
sample quality control, differential expression analysis using
**DESeq2**, **edgeR**, or **limma**, Gene Ontology enrichment,
publication-ready visualization, and automated HTML/PDF/PPTX reporting
within a single reproducible workflow.

The package is designed for both bioinformaticians and experimental
biologists, minimizing manual scripting while preserving analytical
flexibility.

## Features

- Automated bulk RNA-seq workflow
- Raw and post-filtering quality control
- Differential expression with DESeq2, edgeR, or limma
- Pairwise and single-comparison analyses
- Gene Ontology enrichment
- Publication-ready PCA, volcano plots, and heatmaps
- Automated HTML, PDF, and PowerPoint reports
- Reproducibility files and session information
- Human, mouse, rat, and custom OrgDb support

## Workflow

![](reference/figures/DEGgo_workflow.png)

## Installation

``` r


install.packages("remotes")

# Development version
remotes::install_github("ymbouamboua/DEGgo")
```

``` r


library(DEGgo)
```

## Quick start

The package includes a ready-to-use RNA-seq dataset derived from the
Bioconductor **airway** experiment.

``` r


counts <- read.delim(
  system.file("extdata", "airway_counts.tsv", package = "DEGgo"),
  check.names = FALSE
)

metadata <- read.delim(
  system.file("extdata", "airway_metadata.tsv", package = "DEGgo"),
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
  contrast = c("treatment", "PAMH", "PBS"),
  output_dir = "DEGgo_results"
)

results$summary
```

## Pairwise differential expression analysis

``` r


pairwise_contrasts <- list(
  TreatmentA_vs_Control_Male = c(
    "comparison_group",
    "TreatmentA_Male",
    "Control_Male"
  ),
  TreatmentA_vs_Control_Female = c(
    "comparison_group",
    "TreatmentA_Female",
    "Control_Female"
  ),
  Male_vs_Female_Control = c(
    "comparison_group",
    "Control_Male",
    "Control_Female"
  )
)

results <- run_deggo(
  counts = counts,
  metadata = metadata,
  organism = "mouse",
  analysis_mode = "pairwise",
  sample_col = "sample",
  pairwise_group_cols = c("condition", "sex"),
  pairwise_contrasts = pairwise_contrasts,
  output_dir = "DEGgo_results"
)
```

## Functions

| Function | Description |
|:---|:---|
| [`check_raw_counts()`](https://ymbouamboua.github.io/DEGgo/reference/check_raw_counts.md) | Validate raw count matrices and sample identifiers before analysis |
| [`explore_bulk_rnaseq()`](https://ymbouamboua.github.io/DEGgo/reference/explore_bulk_rnaseq.md) | Perform exploratory quality control on raw RNA-seq data |
| [`remove_flagged_samples()`](https://ymbouamboua.github.io/DEGgo/reference/remove_flagged_samples.md) | Remove low-quality samples identified during quality control |
| [`marker_score_check()`](https://ymbouamboua.github.io/DEGgo/reference/marker_score_check.md) | Assess marker gene expression for biological sample validation |
| [`plot_gene_heatmap()`](https://ymbouamboua.github.io/DEGgo/reference/plot_gene_heatmap.md) | Visualize marker or selected genes as a heatmap |
| [`run_deggo()`](https://ymbouamboua.github.io/DEGgo/reference/run_deggo.md) | Run the complete bulk RNA-seq differential expression workflow |
| [`run_sample_qc()`](https://ymbouamboua.github.io/DEGgo/reference/run_sample_qc.md) | Perform post-filtering sample quality control and diagnostics |
| [`run_go_enrichment()`](https://ymbouamboua.github.io/DEGgo/reference/run_go_enrichment.md) | Perform Gene Ontology enrichment analysis |
| [`plot_go_terms()`](https://ymbouamboua.github.io/DEGgo/reference/plot_go_terms.md) | Generate publication-ready Gene Ontology enrichment plots |
| [`plot_all_go_terms()`](https://ymbouamboua.github.io/DEGgo/reference/plot_all_go_terms.md) | Compare Gene Ontology enrichment across multiple analyses |
| [`extract_expression()`](https://ymbouamboua.github.io/DEGgo/reference/extract_expression.md) | Extract raw, normalized, log2-normalized, or VST expression values |
| [`plot_gene_expression()`](https://ymbouamboua.github.io/DEGgo/reference/plot_gene_expression.md) | Visualize expression of selected genes across samples |
| [`plot_heatmap()`](https://ymbouamboua.github.io/DEGgo/reference/plot_heatmap.md) | Generate heatmaps of differentially expressed genes |
| [`plot_pca()`](https://ymbouamboua.github.io/DEGgo/reference/plot_pca.md) | Perform principal component analysis (PCA) |
| [`plot_volcano()`](https://ymbouamboua.github.io/DEGgo/reference/plot_volcano.md) | Generate publication-ready volcano plots |
| [`generate_deggo_report()`](https://ymbouamboua.github.io/DEGgo/reference/generate_deggo_report.md) | Generate an automated HTML or PDF analysis report |
| [`generate_deggo_pptx()`](https://ymbouamboua.github.io/DEGgo/reference/generate_deggo_pptx.md) | Generate a PowerPoint summary report |
| [`deggo_extract_deg_genes()`](https://ymbouamboua.github.io/DEGgo/reference/deggo_extract_deg_genes.md) | Extract selected genes from differential expression results |
| [`deggo_extract_go_genes()`](https://ymbouamboua.github.io/DEGgo/reference/deggo_extract_go_genes.md) | Extract genes associated with enriched Gene Ontology terms |
| [`deggo_extract_go_genes_pairwise()`](https://ymbouamboua.github.io/DEGgo/reference/deggo_extract_go_genes_pairwise.md) | Extract GO-associated genes across pairwise comparisons |
| [`deggo_extract_go_keywords()`](https://ymbouamboua.github.io/DEGgo/reference/deggo_extract_go_keywords.md) | Search enriched Gene Ontology terms using biological keywords |

## Supported organisms

DEGgo provides built-in annotation support for the following organisms:

| Organism                  | Parameter | OrgDb package  |
|:--------------------------|:----------|:---------------|
| Human (*Homo sapiens*)    | `"human"` | `org.Hs.eg.db` |
| Mouse (*Mus musculus*)    | `"mouse"` | `org.Mm.eg.db` |
| Rat (*Rattus norvegicus*) | `"rat"`   | `org.Rn.eg.db` |

Additional organisms can be analyzed by supplying a compatible
Bioconductor `OrgDb` annotation database.

## Documentation

A complete tutorial and advanced examples are available in the package
vignette.

``` r


browseVignettes("DEGgo")
```

## Citation

If you use **DEGgo** in your research, please cite:

> Yvon Mbouamboua, Vincent Prevot, Paolo Giacobini (2026). **DEGgo:
> automated bulk RNA-seq differential expression analysis and Gene
> Ontology enrichment.** Zenodo.
> <https://doi.org/10.5281/zenodo.20785178>

# License

MIT © Yvon MBOUAMBOUA
