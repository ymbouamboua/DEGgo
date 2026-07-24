# DEGgo

![DEGgo logo](../reference/figures/DEGgo_logo.svg)

  
  

**An integrated R framework for automated bulk RNA-seq differential
expression, functional enrichment, circadian rhythmicity analysis, and
reproducible reporting.**

### Overview

**DEGgo** is an open-source R package for reproducible bulk RNA-seq
analysis.

It integrates:

- input validation and sample matching;
- raw and post-filtering quality control;
- differential expression analysis with DESeq2, edgeR, or limma;
- single and pairwise contrasts;
- PCA, volcano plots, heatmaps, and gene-expression visualization;
- Gene Ontology enrichment;
- circadian rhythmicity analysis using MetaCycle and cosinor regression;
- automated HTML, PDF, and PowerPoint reporting;
- reproducibility files, analysis manifests, and session information.

### Installation

``` r

install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github("ymbouamboua/DEGgo")
```

``` r

library(DEGgo)
```

Optional rhythmicity dependencies:

``` r

install.packages(c("MetaCycle", "cosinor", "cosinor2"))
```

### Workflow

``` text
Raw counts
      │
      ▼
 Input validation
      │
      ▼
 Sample QC
      │
      ▼
 Expression filtering
      │
      ▼
 Differential expression
      │
      ├── PCA
      ├── Volcano plots
      ├── Heatmaps
      ├── GO enrichment
      └── Reports
      │
      ▼
 Optional rhythmicity analysis
      ├── MetaCycle
      ├── Cosinor
      ├── Differential rhythmicity
      └── Publication figures
```

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
is specified with `sample_col`.

``` text
sample      condition    batch
Sample1     control      A
Sample2     treated      A
Sample3     control      B
```

Sample identifiers in the metadata must match the sample columns in the
count table.

## Differential expression workflow

### Example data

DEGgo includes a ready-to-use RNA-seq dataset derived from the
Bioconductor [airway](https://bioconductor.org/packages/airway/)
package.

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

### Single comparison mode

``` r

metadata$condition <- metadata$dex

results <- run_deggo(
  counts = counts,
  metadata = metadata,
  project_name = "Airway dexamethasone analysis",
  gene_col = "gene_id",
  organism = "human",
  sample_col = "SampleName",
  method = "DESeq2",
  analysis_mode = "single",
  design_formula = ~ cell + condition,
  contrast = c("condition", "trt", "untrt"),
  filter_method = "count",
  min_count = 5,
  min_samples = 2,
  min_total = 10,
  padj_cutoff = 0.05,
  logfc_cutoff = 0.25,
  output_dir = "DEGgo_airway",
  generate_report = TRUE,
  report_formats = "html",
  generate_pptx = TRUE
)
```

Positive `log2FoldChange` values represent higher expression in the
treated group.

### Pairwise comparison mode

``` r

metadata_mouse$group <- interaction(
  metadata_mouse$treatment,
  metadata_mouse$sex,
  metadata_mouse$tissue,
  sep = "_",
  drop = TRUE
)
```

``` r

pairwise_contrasts <- list(
  treatment_vs_control_male = c(
    "group",
    "treatment_male",
    "control_male"
  ),

  treatment_vs_control_female = c(
    "group",
    "treatment_female",
    "control_female"
  ),

  male_vs_female_control = c(
    "group",
    "control_male",
    "control_female"
  )
)
```

``` r

pairwise_results <- run_deggo(
  counts = counts_mouse,
  metadata = metadata_mouse,
  project_name = "Pairwise analysis",
  organism = "mouse",
  method = "DESeq2",
  analysis_mode = "pairwise",
  sample_col = "sample",
  design_formula = ~ group,
  pairwise_group_cols = c("treatment", "sex", "tissue"),
  pairwise_contrast_col = "group",
  pairwise_contrasts = pairwise_contrasts,
  output_dir = "DEGgo_pairwise_results",
  generate_report = TRUE,
  report_formats = "html",
  generate_pptx = TRUE
)
```

## Quality control

``` r

qc <- explore_bulk_rnaseq(
  counts = counts,
  metadata = metadata,
  sample_col = "SampleName",
  markers = c("TSPAN6", "DPM1", "SCYL3"),
  output_dir = "DEGgo_airway_QC_raw"
)
```

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
```

## Gene expression extraction

``` r

expr_df <- extract_expression(
  dds = results$dds,
  metadata = results$metadata,
  genes = c("TSPAN6", "DPM1", "SCYL3"),
  assay = "vst",
  gene_col = "SYMBOL"
)
```

``` r

plot_gene_expression(
  expr_df,
  gene = "TSPAN6",
  x = "condition",
  color = "condition",
  facet = "cell",
  geom = "violin"
)
```

## Gene Ontology enrichment

``` r

go <- run_go_enrichment(
  sig_deg = results$sig_deg,
  comparison = "condition_trt_vs_untrt",
  ontology = "BP",
  orgdb = org.Hs.eg.db::org.Hs.eg.db,
  output_dir = file.path("DEGgo_airway", "GO_condition_trt_vs_untrt")
)
```

``` r

plot_go_terms(
  go_df = go$go_results,
  comparison = "Dexamethasone treated vs untreated",
  top_n = 15
)
```

## Circadian rhythmicity analysis

### Public circadian example

DEGgo includes a reproducible workflow based on GSE98965 from Mure *et
al.* (Science, 2018).

``` r

public_rhythm <- run_public_circadian_example(
  tissue_code = "LIV",
  output_dir = "Circadian_results"
)
```

``` r

gse98965 <- download_gse98965(
  output_dir = "data/GSE98965"
)

liver_data <- prepare_circadian_tissue(
  expression_data = gse98965,
  tissue_code = "LIV"
)
```

### Run rhythmicity analysis

``` r

rhythm_results <- run_deggo_rhythmicity(
  expr = expression_matrix,
  metadata = metadata,
  sample_col = "sample",
  time_col = "time",
  assay = "log2_normalized",
  methods = c("meta2d", "cosinor"),
  period_range = c(20, 28),
  cycle_length = 24,
  padj_cutoff = 0.05,
  output_dir = "DEGgo_rhythmicity_results"
)
```

### Extract rhythmic genes

``` r

rhythmic_genes <- deggo_extract_rhythmic_genes(
  rhythm_results,
  rhythmic_by = "both",
  padj_cutoff = 0.05
)
```

### Differential rhythmicity

``` r

group_rhythm <- run_deggo_rhythmicity(
  expr = expression_matrix,
  metadata = metadata,
  sample_col = "sample",
  time_col = "time",
  group_col = "condition",
  assay = "log2_normalized",
  methods = c("meta2d", "cosinor"),
  cycle_length = 24,
  padj_cutoff = 0.05,
  output_dir = "DEGgo_group_rhythmicity"
)
```

## Main functions

| Function | Description |
|:---|:---|
| [`run_deggo()`](https://ymbouamboua.github.io/DEGgo/reference/run_deggo.md) | Run the complete differential expression workflow |
| [`run_deggo_rhythmicity()`](https://ymbouamboua.github.io/DEGgo/reference/run_deggo_rhythmicity.md) | Detect rhythmic genes with MetaCycle and cosinor regression |
| [`deggo_extract_rhythmic_genes()`](https://ymbouamboua.github.io/DEGgo/reference/deggo_extract_rhythmic_genes.md) | Extract rhythmic genes from rhythmicity results |
| [`download_gse98965()`](https://ymbouamboua.github.io/DEGgo/reference/download_gse98965.md) | Download the public GSE98965 circadian dataset |
| [`prepare_circadian_tissue()`](https://ymbouamboua.github.io/DEGgo/reference/prepare_circadian_tissue.md) | Prepare one tissue from the public circadian atlas |
| [`run_public_circadian_example()`](https://ymbouamboua.github.io/DEGgo/reference/run_public_circadian_example.md) | Run the complete public circadian example |
| [`check_raw_counts()`](https://ymbouamboua.github.io/DEGgo/reference/check_raw_counts.md) | Validate counts and sample identifiers |
| [`explore_bulk_rnaseq()`](https://ymbouamboua.github.io/DEGgo/reference/explore_bulk_rnaseq.md) | Perform exploratory quality control |
| [`remove_flagged_samples()`](https://ymbouamboua.github.io/DEGgo/reference/remove_flagged_samples.md) | Remove low-quality samples |
| [`run_sample_qc()`](https://ymbouamboua.github.io/DEGgo/reference/run_sample_qc.md) | Perform post-filtering sample quality control |
| [`run_go_enrichment()`](https://ymbouamboua.github.io/DEGgo/reference/run_go_enrichment.md) | Perform Gene Ontology enrichment |
| [`plot_heatmap()`](https://ymbouamboua.github.io/DEGgo/reference/plot_heatmap.md) | Plot differentially expressed genes |
| [`plot_pca()`](https://ymbouamboua.github.io/DEGgo/reference/plot_pca.md) | Generate a PCA plot |
| [`plot_volcano()`](https://ymbouamboua.github.io/DEGgo/reference/plot_volcano.md) | Generate a volcano plot |
| [`generate_deggo_report()`](https://ymbouamboua.github.io/DEGgo/reference/generate_deggo_report.md) | Generate HTML or PDF reports |
| [`generate_deggo_pptx()`](https://ymbouamboua.github.io/DEGgo/reference/generate_deggo_pptx.md) | Generate a PowerPoint summary |
| [`deggo_palette()`](https://ymbouamboua.github.io/DEGgo/reference/deggo_palette.md) | Retrieve DEGgo colour palettes |

## Supported organisms

| Organism                  | Parameter | OrgDb package  |
|:--------------------------|:----------|:---------------|
| Human (*Homo sapiens*)    | `"human"` | `org.Hs.eg.db` |
| Mouse (*Mus musculus*)    | `"mouse"` | `org.Mm.eg.db` |
| Rat (*Rattus norvegicus*) | `"rat"`   | `org.Rn.eg.db` |

Additional organisms can be analyzed with a compatible Bioconductor
`OrgDb` object.

## Citation

If you use DEGgo, please cite the Zenodo release:

<https://doi.org/10.5281/zenodo.20785178>

## Contributing

Bug reports, feature requests, and pull requests are welcome through the
[DEGgo GitHub repository](https://github.com/ymbouamboua/DEGgo).

## License

MIT License © Yvon MBOUAMBOUA
