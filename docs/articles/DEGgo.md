# DEGgo

![DEGgo package logo showing an integrated workflow for bulk RNA-seq
analysis](../reference/figures/DEGgo_logo.png)

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

The package can be installed from GitHub:

``` r

install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github("ymbouamboua/DEGgo")
```

Load DEGgo:

``` r

library(DEGgo)
```

### Workflow

![DEGgo bulk RNA-seq workflow from raw counts and metadata through
quality control, differential expression analysis, Gene Ontology
enrichment, visualization and HTML report
generation](../reference/figures/DEGgo_workflow.png)

### Input data

#### Count table

DEGgo accepts raw count tables or matrices. Gene annotation columns can
be automatically detected using `gene_col` and `feature_col`.

``` text
gene_id      gene_name    Sample1    Sample2    Sample3
ENSMUSG1     Adipoq       120        145        98
ENSMUSG2     Lep          65         80         50
ENSMUSG3     Ucp1         12         18         250
```

#### Metadata

Metadata must contain one row per sample. A `sample` column is expected
or automatically detected using `sample_col`.

``` text
sample      tissue    sex       treatment
Sample1     WAT       Female    PBS
Sample2     WAT       Female    PAMH
Sample3     BAT       Female    PBS
```

For pairwise analysis, DEGgo can build a comparison variable from
multiple metadata columns:

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

### Main functions

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

## Complete DEGgo workflow

### 1. Raw count quality control

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

### 2. Remove flagged samples

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

### 3. Re-run QC after cleaning

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

### 4. Tissue marker validation

``` r

marker_sets <- list(
  BAT = c("Ucp1", "Cidea", "Ppargc1a", "Prdm16", "Dio2", "Elovl3"),
  WAT = c("Adipoq", "Lep", "Retn", "Fabp4", "Pparg", "Lpl"),
  TESTIS = c("Amh", "Amhr2", "Sox9", "Dhh", "Star", "Cyp11a1", "Cyp17a1", "Insl3"),
  OVARY = c("Foxl2", "Fshr", "Cyp19a1", "Bmp15", "Gdf9")
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

## Differential expression analysis

### Single comparison mode

Use `analysis_mode = "single"` when you have one contrast such as `PAMH`
vs `PBS`, `KO` vs `WT`, or `treated` vs `control`.

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

### Pairwise contrast mode

Use `analysis_mode = "pairwise"` when your experiment contains multiple
structured comparisons, such as treatment by sex by tissue.

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

## Automated HTML/PDF report

DEGgo can automatically generate an HTML or PDF report at the end of
[`run_deggo()`](https://ymbouamboua.github.io/DEGgo/reference/run_deggo.md).

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

## GO enrichment

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

``` r

plot_go_terms(
  go_df = go$go_results,
  comparison = "WAT Female PAMH vs PBS",
  top_n = 15
)
```

## Extract genes and GO terms of interest

``` r

genes_interest <- c(
  "Amh", "Amhr2", "Ar", "Akr1c3",
  "Star", "Cyp11a1", "Cyp17a1", "Cyp19a1",
  "Foxl2", "Fshr", "Sox9", "Insl3",
  "Pparg", "Fabp4", "Adipoq", "Lep",
  "Ucp1", "Prdm16", "Ppargc1a", "Cidea"
)
```

``` r

gene_summary <- deggo_extract_deg_genes(results, genes_interest)
gene_summary
```

``` r

go_keywords <- c(
  "tgf", "bmp", "smad", "inflammation", "immune",
  "steroidogenesis", "androgen", "steroid",
  "adipocyte", "lipid", "fatty acid",
  "thermogenesis", "mitochond", "hormone",
  "reproductive", "gonad", "ovary", "testis"
)

go_keywords_results <- deggo_extract_go_keywords(
  results,
  keywords = go_keywords
)

go_keywords_results
```

## Output structure

``` text
DEGgo_results/
└── DEGgo_YYYY_MM_DD_DESeq2_pairwise/
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

results <- run_deggo(
  counts = counts,
  metadata = metadata,
  organism = "custom",
  orgdb = org.Custom.eg.db
)
```

## Platform compatibility and HPC support

DEGgo is implemented entirely in R and is compatible with Linux, macOS,
and Windows. It can be used on local computers, HPC clusters, and
containerized environments such as Docker or Apptainer/Singularity.

## Citation

If you use DEGgo in your work, please cite:

> MBOUAMBOUA Y., Prevot V., Giacobini P. DEGgo: an R package for
> automated bulk RNA-seq differential expression analysis and Gene
> Ontology enrichment.

## License

MIT License © Yvon MBOUAMBOUA
