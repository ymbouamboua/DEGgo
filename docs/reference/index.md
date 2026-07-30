# Package index

## Main workflow

Core functions for running differential expression analysis, functional
enrichment, report generation, and result export.

- [`run_deggo()`](https://ymbouamboua.github.io/DEGgo/reference/run_deggo.md)
  : Run DEGgo bulk RNA-seq downstream analysis
- [`run_go_enrichment()`](https://ymbouamboua.github.io/DEGgo/reference/run_go_enrichment.md)
  : Run Gene Ontology enrichment analysis
- [`generate_deggo_report()`](https://ymbouamboua.github.io/DEGgo/reference/generate_deggo_report.md)
  : Generate a DEGgo report
- [`generate_deggo_pptx()`](https://ymbouamboua.github.io/DEGgo/reference/generate_deggo_pptx.md)
  : Generate a PowerPoint report from DEGgo results

## Experimental design assessment

Functions for validating RNA-seq study designs, identifying repeated
experimental units, assessing dominant sources of global expression
variation, detecting potential pseudoreplication, and recommending
appropriate DESeq2 or dream models.

- [`design_qc()`](https://ymbouamboua.github.io/DEGgo/reference/design_qc.md)
  : Experimental design quality control for bulk RNA-seq
- [`aggregate_counts_by_unit()`](https://ymbouamboua.github.io/DEGgo/reference/aggregate_counts_by_unit.md)
  : Aggregate raw counts by an experimental unit
- [`plot(`*`<deggo_design_qc>`*`)`](https://ymbouamboua.github.io/DEGgo/reference/plot.deggo_design_qc.md)
  : Plot DesignQC PCA
- [`print(`*`<deggo_design_qc>`*`)`](https://ymbouamboua.github.io/DEGgo/reference/print.deggo_design_qc.md)
  : Print a DEGgo design quality control summary

## Circadian rhythmicity analysis

Functions for preparing circadian transcriptomic datasets, detecting
rhythmic genes, testing differential rhythmicity, and running
reproducible public examples.

- [`run_deggo_rhythmicity()`](https://ymbouamboua.github.io/DEGgo/reference/run_deggo_rhythmicity.md)
  : Run DEGgo rhythmicity (circadian) analysis
- [`deggo_extract_rhythmic_genes()`](https://ymbouamboua.github.io/DEGgo/reference/deggo_extract_rhythmic_genes.md)
  : Extract selected genes from DEGgo rhythmicity results
- [`prepare_circadian_tissue()`](https://ymbouamboua.github.io/DEGgo/reference/prepare_circadian_tissue.md)
  : Prepare one GSE98965 tissue for circadian analysis
- [`download_gse98965()`](https://ymbouamboua.github.io/DEGgo/reference/download_gse98965.md)
  : Download the GSE98965 baboon transcriptomic atlas
- [`run_public_circadian_example()`](https://ymbouamboua.github.io/DEGgo/reference/run_public_circadian_example.md)
  : Run the public GSE98965 circadian RNA-seq example

## Quality control

Functions for count-data exploration, sample-level quality control,
filtering, sample removal, and marker assessment.

- [`explore_bulk_rnaseq()`](https://ymbouamboua.github.io/DEGgo/reference/explore_bulk_rnaseq.md)
  : Explore and QC bulk RNA-seq count data
- [`check_raw_counts()`](https://ymbouamboua.github.io/DEGgo/reference/check_raw_counts.md)
  : Validate raw count table
- [`run_sample_qc()`](https://ymbouamboua.github.io/DEGgo/reference/run_sample_qc.md)
  : Run sample-level RNA-seq quality control
- [`remove_flagged_samples()`](https://ymbouamboua.github.io/DEGgo/reference/remove_flagged_samples.md)
  : Remove flagged samples from counts and metadata
- [`marker_score_check()`](https://ymbouamboua.github.io/DEGgo/reference/marker_score_check.md)
  : Validate sample identity using marker gene signatures

## Visualization

Functions for generating publication-ready differential expression,
gene-expression, and functional enrichment visualizations.

- [`plot_volcano()`](https://ymbouamboua.github.io/DEGgo/reference/plot_volcano.md)
  : Generate volcano plot for differential expression results
- [`plot_heatmap()`](https://ymbouamboua.github.io/DEGgo/reference/plot_heatmap.md)
  : Plot a DEGgo Differential Expression Heatmap
- [`plot_pca()`](https://ymbouamboua.github.io/DEGgo/reference/plot_pca.md)
  : Generate DEGgo PCA plot
- [`plot_go_terms()`](https://ymbouamboua.github.io/DEGgo/reference/plot_go_terms.md)
  : Plot GO terms by regulation status
- [`plot_all_go_terms()`](https://ymbouamboua.github.io/DEGgo/reference/plot_all_go_terms.md)
  : Plot GO enrichment terms across all DEGgo comparisons
- [`plot_gene_expression()`](https://ymbouamboua.github.io/DEGgo/reference/plot_gene_expression.md)
  : Plot normalized gene expression
- [`plot_gene_heatmap()`](https://ymbouamboua.github.io/DEGgo/reference/plot_gene_heatmap.md)
  : Plot Expression Heatmap for Selected Genes

## Expression analysis

Functions for extracting and inspecting normalized or transformed
gene-expression values.

- [`extract_expression()`](https://ymbouamboua.github.io/DEGgo/reference/extract_expression.md)
  : Extract normalized gene expression

## Gene and GO extraction

Functions for extracting differentially expressed genes, enriched Gene
Ontology terms, rhythmic genes, and associated gene sets.

- [`deggo_extract_deg_genes()`](https://ymbouamboua.github.io/DEGgo/reference/deggo_extract_deg_genes.md)
  : Extract selected DEG genes from DEGgo results
- [`deggo_extract_go_genes()`](https://ymbouamboua.github.io/DEGgo/reference/deggo_extract_go_genes.md)
  : Extract genes found in GO terms
- [`deggo_extract_go_genes_pairwise()`](https://ymbouamboua.github.io/DEGgo/reference/deggo_extract_go_genes_pairwise.md)
  : Extract selected genes from pairwise GO results
- [`deggo_extract_go_keywords()`](https://ymbouamboua.github.io/DEGgo/reference/deggo_extract_go_keywords.md)
  : Extract GO terms matching keywords

## Palettes and visualization helpers

Colour palettes and visual styling utilities used throughout DEGgo.

- [`deggo_palette()`](https://ymbouamboua.github.io/DEGgo/reference/deggo_palette.md)
  : Retrieve a DEGgo color palette
