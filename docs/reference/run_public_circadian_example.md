# Run the public GSE98965 circadian RNA-seq example

Run a complete, reproducible circadian rhythmicity analysis using one
tissue from the public GSE98965 baboon transcriptomic atlas.

## Usage

``` r
run_public_circadian_example(
  tissue_code = "LIV",
  data_dir = "GSE98965",
  output_dir = "DEGgo_rhythmicity_results",
  assay = c("log2_normalized", "vst", "raw", "normalized"),
  min_expression = 1,
  min_samples = 3,
  padj_cutoff = 0.05,
  n_top_plots = 25,
  show_gene_id = TRUE,
  seed = 4173
)
```

## Arguments

- tissue_code:

  Character string specifying the GSE98965 tissue code. The default
  `"LIV"` corresponds to liver.

- data_dir:

  Character string. Directory used to store the downloaded GSE98965
  supplementary files.

- output_dir:

  Character string. Parent output directory for the analysis. A
  tissue-specific subdirectory is created automatically.

- assay:

  Character string indicating the expression representation supplied to
  [`run_deggo_rhythmicity()`](https://ymbouamboua.github.io/DEGgo/reference/run_deggo_rhythmicity.md).
  One of `"log2_normalized"`, `"vst"`, `"raw"`, or `"normalized"`.

  The default is `"log2_normalized"` because
  [`prepare_circadian_tissue()`](https://ymbouamboua.github.io/DEGgo/reference/prepare_circadian_tissue.md)
  returns `log2(FPKM + 1)` values.

- min_expression:

  Numeric value. Minimum FPKM expression used to define an expressed
  gene.

- min_samples:

  Positive integer. Minimum number of tissue samples in which a gene
  must reach `min_expression`.

- padj_cutoff:

  Numeric value between zero and one. Adjusted p-value threshold used to
  classify rhythmic genes.

- n_top_plots:

  Positive integer. Number of top rhythmic genes for which individual
  expression and fitted-rhythm plots are generated.

- show_gene_id:

  Logical. If `TRUE`, rhythmicity plot labels include both the gene
  symbol and matrix gene identifier.

- seed:

  Integer random seed used for reproducibility.

## Value

A named list containing:

- summary:

  Combined rhythmicity result table.

- rhythmic_genes:

  Genes classified as rhythmic by at least one selected method.

- rhythmicity_summary:

  Counts of genes by rhythmicity class.

- analysis_summary:

  Summary of the public dataset analysis.

- prepared:

  Prepared expression matrix, metadata, and gene annotation.

- deggo:

  Complete object returned by
  [`run_deggo_rhythmicity()`](https://ymbouamboua.github.io/DEGgo/reference/run_deggo_rhythmicity.md).

- output_dir:

  Tissue-specific output directory.

## Details

The workflow downloads the public processed expression matrix, selects
one tissue, filters lowly expressed genes, converts Zeitgeber Time
labels into numeric values, validates expression-metadata
correspondence, runs MetaCycle and cosinor analysis, adds gene-symbol
annotations, generates diagnostic plots, and exports result tables.

The default analysis uses:

- MetaCycle algorithms ARS, JTK, and Lomb-Scargle;

- a tested period range of 20 to 28 hours;

- a 24-hour single-component cosinor model;

- Benjamini-Hochberg adjusted p-values;

- gene symbols for output tables and figure labels.

The workflow is intended as a reproducible demonstration and as a
starting point for custom circadian RNA-seq analyses.

## References

Mure LS et al. Diurnal transcriptome atlas of a primate across major
neural and peripheral tissues. Science. 2018.

GEO accession: GSE98965.

## See also

[`download_gse98965()`](https://ymbouamboua.github.io/DEGgo/reference/download_gse98965.md),
[`prepare_circadian_tissue()`](https://ymbouamboua.github.io/DEGgo/reference/prepare_circadian_tissue.md),
[`run_deggo_rhythmicity()`](https://ymbouamboua.github.io/DEGgo/reference/run_deggo_rhythmicity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
results <- run_public_circadian_example(
  tissue_code = "LIV",
  data_dir = "GSE98965",
  output_dir = "DEGgo_rhythmicity_results",
  min_expression = 1,
  min_samples = 3,
  padj_cutoff = 0.05,
  n_top_plots = 25,
  show_gene_id = FALSE,
  seed = 4173
)

head(results$summary)
head(results$rhythmic_genes)
results$rhythmicity_summary
} # }
```
