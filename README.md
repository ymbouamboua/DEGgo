DEGgo
================

- [Overview](#overview)
- [Workflow](#workflow)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Main functions](#main-functions)
- [Supported organisms](#supported-organisms)
- [Documentation](#documentation)
- [Citation](#citation)
- [License](#license)

<p align="center">

<img src="https://img.shields.io/badge/R-%3E%3D4.3-blue">
<img src="https://img.shields.io/badge/license-MIT-green">
<img src="https://img.shields.io/github/stars/ymbouamboua/DEGgo?style=social">
</p>

<p align="center">

<img src="man/figures/DEGgo_logo.png"
       alt="DEGgo package logo showing an integrated workflow for bulk RNA-seq analysis"
       width="350">
</p>

## Overview

DEGgo is an R package for automated bulk RNA-seq downstream analysis.

It provides an end-to-end workflow from raw count matrices and sample
metadata to:

- quality control;
- sample validation;
- differential expression analysis;
- Gene Ontology enrichment;
- publication-ready visualizations;
- automated HTML/PDF reporting.

DEGgo supports DESeq2, edgeR, and limma and is designed for both
bioinformaticians and experimental biologists.

## Workflow

<p align="center">

<img src="man/figures/DEGgo_workflow.png" width="600">
</p>

## Installation

``` r

install.packages("remotes")
remotes::install_github("ymbouamboua/DEGgo")
```

``` r
library(DEGgo)
```

## Quick start

``` r

results <- run_deggo(
  counts = counts,
  metadata = metadata,
  organism = "mouse",
  method = "DESeq2",
  analysis_mode = "single",
  contrast = c("condition", "treated", "control")
)

results$summary
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

## Supported organisms

DEGgo provides built-in annotation support for the following organisms:

| Organism                  | Parameter | OrgDb package  |
|:--------------------------|:----------|:---------------|
| Human (*Homo sapiens*)    | `"human"` | `org.Hs.eg.db` |
| Mouse (*Mus musculus*)    | `"mouse"` | `org.Mm.eg.db` |
| Rat (*Rattus norvegicus*) | `"rat"`   | `org.Rn.eg.db` |

Custom organisms are supported through user-supplied Bioconductor OrgDb
annotation databases.

## Documentation

The complete DEGgo tutorial and advanced workflows are available in the
package vignette:

``` r

browseVignettes("DEGgo")
```

# Citation

If you use DEGgo in your work, please cite:

> Yvon MBOUAMBOUA, Vincent Prevot and Paolo Giacobini. DEGgo: an
> integrated framework for bulk RNA-seq differential expression analysis
> and functional enrichment.

# License

MIT © Yvon MBOUAMBOUA
