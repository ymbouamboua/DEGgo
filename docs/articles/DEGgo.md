# DEGgo

![DEGgo package logo](../reference/figures/DEGgo_logo.png)

![R version badge](https://img.shields.io/badge/R-%3E%3D4.3-blue)![MIT
license badge](https://img.shields.io/badge/license-MIT-green)![GitHub
stars
badge](https://img.shields.io/github/stars/ymbouamboua/DEGgo?style=social)

### Overview

**DEGgo** is an R package for end-to-end bulk RNA-seq differential
expression analysis and Gene Ontology enrichment.

DEGgo provides a complete workflow from raw count matrices and sample
metadata to biological interpretation:

- raw count quality control;
- sample validation and filtering;
- tissue marker verification;
- differential expression analysis;
- single and pairwise contrasts;
- Gene Ontology enrichment;
- extraction of genes and GO terms of interest;
- publication-ready visualization;
- automated HTML/PDF reporting.

### Installation

Install DEGgo from GitHub:

``` r

install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github("ymbouamboua/DEGgo")
```

Load the package:

``` r

library(DEGgo)
```

### Workflow

![DEGgo bulk RNA-seq workflow](../reference/figures/DEGgo_workflow.png)

DEGgo follows a standard bulk RNA-seq workflow:

1.  input preparation and sample matching;
2.  quality control;
3.  low-expression filtering;
4.  differential expression analysis;
5.  DEG annotation and export;
6.  PCA, heatmap and volcano visualization;
7.  Gene Ontology enrichment;
8.  automated report generation.

### Input data

#### Count table

DEGgo accepts raw count tables or matrices. The count table should
contain one gene identifier column and one column per sample.

``` text
gene_id          gene_name    Sample1    Sample2    Sample3
ENSG00000000003  TSPAN6       120        145        98
ENSG00000000005  TNMD         65         80         50
ENSG00000000419  DPM1         12         18         250
```

#### Metadata

Metadata must contain one row per sample. The sample identifier column
is supplied using `sample_col`.

``` text
sample      condition    batch
Sample1     control      A
Sample2     treated      A
Sample3     control      B
```

The sample names in the metadata must match the sample columns in the
count table.

## Example analysis with the airway dataset

This vignette demonstrates a complete DEGgo workflow using a small
DEGgo-ready dataset derived from the Bioconductor `airway` package.

The dataset contains human RNA-seq raw counts from airway smooth muscle
cells treated with dexamethasone. The differential expression analysis
compares treated samples (`trt`) against untreated samples (`untrt`)
while accounting for the cell line effect.

### Load example data

The example data are stored in `inst/extdata` and can be accessed with
[`system.file()`](https://rdrr.io/r/base/system.file.html).

``` r

counts <- read.delim(
  system.file("extdata", "airway_counts.tsv", package = "DEGgo"),
  check.names = FALSE
)

metadata <- read.delim(
  system.file("extdata", "airway_metadata.tsv", package = "DEGgo"),
  check.names = FALSE
)

dim(counts)
head(counts[, 1:5])
head(metadata)
```

### Run DEGgo on airway

``` r

results <- run_deggo(
  counts = counts,
  metadata = metadata,
  gene_col = "gene_id",
  organism = "human",
  sample_col = "SampleName",
  method = "DESeq2",
  analysis_mode = "single",
  design_formula = ~ cell + dex,
  contrast = c("dex", "trt", "untrt"),
  filter_method = "count",
  min_count = 5,
  min_samples = 2,
  min_total = 10,
  padj_cutoff = 0.05,
  logfc_cutoff = 0.25,
  output_dir = "DEGgo_airway",
  generate_report = TRUE,
  report_formats = "html"
)

results$summary
results$report_files
```

### Output files

DEGgo creates a complete output directory containing:

- differential expression results;
- significant DEG tables;
- volcano plots;
- PCA plots;
- heatmaps;
- Gene Ontology enrichment results;
- an HTML report.

The main report is available at:

``` r

file.path("DEGgo_airway", "DEGgo_Report.html")
```

## Main functions

| Function | Description |
|----|----|
| [`explore_bulk_rnaseq()`](https://ymbouamboua.github.io/DEGgo/reference/explore_bulk_rnaseq.md) | Raw and cleaned RNA-seq quality control |
| [`remove_flagged_samples()`](https://ymbouamboua.github.io/DEGgo/reference/remove_flagged_samples.md) | Remove failed or low-quality samples |
| [`marker_score_check()`](https://ymbouamboua.github.io/DEGgo/reference/marker_score_check.md) | Tissue marker scoring and sample swap detection |
| [`plot_gene_heatmap()`](https://ymbouamboua.github.io/DEGgo/reference/plot_gene_heatmap.md) | Marker or selected-gene heatmap |
| [`run_deggo()`](https://ymbouamboua.github.io/DEGgo/reference/run_deggo.md) | Main DEGgo differential expression workflow |
| [`run_go_enrichment()`](https://ymbouamboua.github.io/DEGgo/reference/run_go_enrichment.md) | Gene Ontology enrichment for DEG results |
| [`plot_go_terms()`](https://ymbouamboua.github.io/DEGgo/reference/plot_go_terms.md) | Publication-ready GO plot |
| [`plot_all_go_terms()`](https://ymbouamboua.github.io/DEGgo/reference/plot_all_go_terms.md) | Plot GO terms across multiple comparisons |
| [`extract_expression()`](https://ymbouamboua.github.io/DEGgo/reference/extract_expression.md) | Extract raw, normalized, log2-normalized or VST expression values |
| [`plot_gene_expression()`](https://ymbouamboua.github.io/DEGgo/reference/plot_gene_expression.md) | Plot normalized gene expression for selected genes |
| [`generate_deggo_report()`](https://ymbouamboua.github.io/DEGgo/reference/generate_deggo_report.md) | Generate HTML/PDF DEGgo report |
| [`deggo_extract_deg_genes()`](https://ymbouamboua.github.io/DEGgo/reference/deggo_extract_deg_genes.md) | Extract genes of interest from DEG results |
| [`deggo_extract_go_genes_pairwise()`](https://ymbouamboua.github.io/DEGgo/reference/deggo_extract_go_genes_pairwise.md) | Extract GO terms containing selected genes |
| [`deggo_extract_go_keywords()`](https://ymbouamboua.github.io/DEGgo/reference/deggo_extract_go_keywords.md) | Extract GO terms matching biological keywords |

## Quality control workflow

Before differential expression analysis, DEGgo can explore library size,
detected genes, sample correlation, clustering and marker expression.

``` r

qc <- explore_bulk_rnaseq(
  counts = counts,
  metadata = metadata,
  markers = c("TSPAN6", "DPM1", "SCYL3"),
  output_dir = "DEGgo_airway_QC_raw"
)

qc$qc
```

Flagged samples can be removed before re-running the analysis.

``` r

cleaned <- remove_flagged_samples(
  counts = counts,
  metadata = metadata,
  qc_table = qc$qc,
  sample_col = "SampleName",
  remove_col = "recommend_remove",
  gene_cols = "gene_id",
  verbose = TRUE
)

counts_clean <- cleaned$counts
metadata_clean <- cleaned$metadata
```

## Single comparison mode

Use `analysis_mode = "single"` when the experiment contains one main
contrast, such as treated vs control.

``` r

results_single <- run_deggo(
  counts = counts,
  metadata = metadata,
  gene_col = "gene_id",
  organism = "human",
  method = "DESeq2",
  analysis_mode = "single",
  sample_col = "SampleName",
  design_formula = ~ cell + dex,
  contrast = c("dex", "trt", "untrt"),
  filter_method = "count",
  min_count = 5,
  min_samples = 2,
  min_total = 10,
  generate_report = TRUE,
  report_formats = "html",
  output_dir = "DEGgo_single"
)

results_single$summary
results_single$report_files
```

## Pairwise contrast mode

Use `analysis_mode = "pairwise"` when the experiment contains multiple
structured comparisons, such as treatment by sex by tissue.

The following example is a template for a mouse experiment with
treatment, sex and tissue metadata columns.

``` r

pairwise_contrasts <- list(
  WAT_Female_PAMH_vs_PBS = c("comparison_group", "PAMH_Female_WAT", "PBS_Female_WAT"),
  WAT_Male_PAMH_vs_PBS = c("comparison_group", "PAMH_Male_WAT", "PBS_Male_WAT"),
  BAT_Female_PAMH_vs_PBS = c("comparison_group", "PAMH_Female_BAT", "PBS_Female_BAT"),
  BAT_Male_PAMH_vs_PBS = c("comparison_group", "PAMH_Male_BAT", "PBS_Male_BAT"),
  Ovary_PAMH_vs_PBS = c("comparison_group", "PAMH_Female_OVARY", "PBS_Female_OVARY"),
  Testis_PAMH_vs_PBS = c("comparison_group", "PAMH_Male_TESTIS", "PBS_Male_TESTIS")
)
```

``` r

results_pairwise <- run_deggo(
  counts = counts_mouse,
  metadata = metadata_mouse,
  organism = "mouse",
  method = "DESeq2",
  analysis_mode = "pairwise",
  sample_col = "sample",
  pairwise_group_cols = c("treatment", "sex", "tissue"),
  pairwise_contrast_col = "comparison_group",
  pairwise_contrasts = pairwise_contrasts,
  filter_method = "count",
  min_count = 2,
  min_samples = 2,
  min_total = 10,
  generate_report = TRUE,
  report_formats = "html",
  output_dir = "DEGgo_pairwise"
)

results_pairwise$summary
results_pairwise$report_files
```

## Automated HTML/PDF report

DEGgo can automatically generate an HTML or PDF report from an existing
result object.

``` r

report_files <- generate_deggo_report(
  results = results,
  output_dir = results$output_dir,
  formats = c("html", "pdf")
)

report_files
```

PDF generation requires a LaTeX installation, for example TinyTeX.

``` r

install.packages("tinytex", repos = "https://cloud.r-project.org")
tinytex::install_tinytex()
```

## Gene expression extraction and plotting

Expression values can be extracted from the fitted object and plotted
for selected genes.

``` r

expr_df <- extract_expression(
  dds = results$dds,
  metadata = results$metadata,
  genes = c("TSPAN6", "DPM1", "SCYL3"),
  assay = "vst",
  gene_col = "SYMBOL"
)

head(expr_df)
```

``` r

plot_gene_expression(
  expr_df,
  gene = "TSPAN6",
  x = "dex",
  color = "dex",
  facet = "cell",
  geom = "violin"
)
```

## GO enrichment

GO enrichment is automatically run inside
[`run_deggo()`](https://ymbouamboua.github.io/DEGgo/reference/run_deggo.md)
when significant genes are detected. It can also be run manually.

``` r

go <- run_go_enrichment(
  sig_deg = results$sig_deg,
  comparison = "dex_trt_vs_untrt",
  ontology = "BP",
  orgdb = org.Hs.eg.db::org.Hs.eg.db,
  output_dir = file.path("DEGgo_airway", "GO_dex_trt_vs_untrt")
)

go$go_results
go$go_plot
```

``` r

plot_go_terms(
  go_df = go$go_results,
  comparison = "Dexamethasone treated vs untreated",
  top_n = 15
)
```

## Extract genes and GO terms of interest

``` r

genes_interest <- c(
  "FKBP5", "DUSP1", "KLF15", "PER1",
  "TSC22D3", "ZBTB16", "TSPAN6", "DPM1"
)
```

``` r

gene_summary <- deggo_extract_deg_genes(results, genes_interest)
gene_summary
```

``` r

go_keywords <- c(
  "glucocorticoid", "steroid", "inflammation",
  "immune", "cytokine", "airway", "epithelial"
)

go_keywords_results <- deggo_extract_go_keywords(
  results,
  keywords = go_keywords
)

go_keywords_results
```

## Output structure

A typical DEGgo output directory contains:

``` text
DEGgo_airway/
├── DEGgo_Report.html
├── DESeq2_results/
├── DESeq2_significant/
├── volcano/
├── heatmaps/
├── PCA/
├── GO/
├── GO_plots/
├── DEGgo_output_manifest.tsv
├── session_info.txt
└── DEGgo_results.rds
```

For pairwise analyses, DEGgo creates comparison-specific directories:

``` text
DEGgo_pairwise/
├── pairwise_DESeq2_results/
├── pairwise_DESeq2_significant/
├── pairwise_volcano/
├── pairwise_heatmaps/
├── pairwise_PCA/
├── pairwise_GO/
├── pairwise_GO_plots/
├── pairwise_summary.tsv
├── DEGgo_Report.html
├── session_info.txt
└── DEGgo_results.rds
```

## Supported organisms

| Organism                  | Parameter | OrgDb package  |
|:--------------------------|:----------|:---------------|
| Human (*Homo sapiens*)    | `"human"` | `org.Hs.eg.db` |
| Mouse (*Mus musculus*)    | `"mouse"` | `org.Mm.eg.db` |
| Rat (*Rattus norvegicus*) | `"rat"`   | `org.Rn.eg.db` |

### Custom annotation databases

Users can analyze additional organisms by supplying a
Bioconductor-compatible OrgDb object through the `orgdb` argument.

``` r

library(org.Custom.eg.db)

results_custom <- run_deggo(
  counts = counts,
  metadata = metadata,
  organism = "custom",
  orgdb = org.Custom.eg.db
)
```

## Platform compatibility and HPC support

DEGgo is implemented entirely in R and is compatible with Linux, macOS
and Windows. It can be used on local computers, HPC clusters and
containerized environments such as Docker or Apptainer/Singularity.

## Citation

If you use DEGgo in your work, please cite:

> MBOUAMBOUA Y., Prevot V., Giacobini P. DEGgo: an R package for
> automated bulk RNA-seq differential expression analysis and Gene
> Ontology enrichment.

## License

MIT License © Yvon MBOUAMBOUA
