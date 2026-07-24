# Run DEGgo rhythmicity (circadian) analysis

Standalone rhythmicity detection workflow using MetaCycle (`meta2d`) and
single-component cosinor regression, with an optional cosinor-based test
for differential rhythmicity between two groups.

## Usage

``` r
run_deggo_rhythmicity(
  expr,
  metadata,
  sample_col = "sample",
  time_col = "time",
  group_col = NULL,
  assay = c("vst", "raw", "normalized", "log2_normalized"),
  methods = c("meta2d", "cosinor"),
  period_range = c(20, 28),
  cycle_length = 24,
  cycMethod = c("ARS", "JTK", "LS"),
  padj_cutoff = 0.05,
  cosinor_engine = c("auto", "package", "manual"),
  output_dir = "DEGgo_rhythmicity",
  project_name = "DEGgo rhythmicity analysis",
  generate_plots = TRUE,
  n_top_plots = 20,
  txtsize = 12,
  seed = 4173,
  gene_annotation = NULL,
  gene_id_col = "gene_id",
  gene_symbol_col = NULL,
  show_gene_id = TRUE,
  verbose = TRUE
)
```

## Arguments

- expr:

  A genes-by-samples numeric matrix/data frame of expression values, or
  a DESeq2 `DESeqDataSet` object.

- metadata:

  Sample metadata data frame. Must contain a sample identifier column
  defined by `sample_col` and a numeric time column defined by
  `time_col`.

- sample_col:

  Column in `metadata` containing sample identifiers that match the
  column names of `expr`.

- time_col:

  Column in `metadata` containing numeric time values, such as hours or
  Zeitgeber Time.

- group_col:

  Optional column in `metadata` defining exactly two groups for
  differential-rhythmicity testing.

- assay:

  Assay to extract when `expr` is a DESeq2 `DESeqDataSet`. One of
  `"vst"`, `"normalized"`, `"log2_normalized"`, or `"raw"`.

- methods:

  Rhythmicity-detection methods to run. One or both of `"meta2d"` and
  `"cosinor"`.

- period_range:

  Numeric vector of length two defining MetaCycle's period search range
  in the same units as `time_col`.

- cycle_length:

  Assumed period for the single-component cosinor model.

- cycMethod:

  Character vector of MetaCycle methods to combine. Allowed values
  include `"ARS"`, `"JTK"`, and `"LS"`.

- padj_cutoff:

  Adjusted p-value cutoff used to classify rhythmic genes.

- cosinor_engine:

  Cosinor fitting engine. `"auto"` uses the `cosinor` and `cosinor2`
  packages when installed and otherwise falls back to the manual
  [`lm()`](https://rdrr.io/r/stats/lm.html) implementation. `"package"`
  requires the optional packages. `"manual"` always uses the internal
  [`lm()`](https://rdrr.io/r/stats/lm.html) implementation.

- output_dir:

  Output directory.

- project_name:

  Optional project name.

- generate_plots:

  Logical. Generate rhythmicity diagnostic plots.

- n_top_plots:

  Number of top rhythmic genes to plot individually.

- txtsize:

  Base text size.

- seed:

  Random seed.

- gene_annotation:

  Optional gene-annotation data frame used to map expression-matrix row
  identifiers to readable gene symbols.

- gene_id_col:

  Gene-identifier column in `gene_annotation`.

- gene_symbol_col:

  Optional gene-symbol column in `gene_annotation`. When `NULL`, common
  symbol-column names are detected automatically.

- show_gene_id:

  Logical. Include the original gene identifier in individual
  rhythmicity plot titles.

- verbose:

  Logical. Print progress messages.

## Value

A `deggo_rhythm_results` object containing the combined summary,
method-specific result tables, matched metadata, expression matrix,
annotation settings, parameters, output paths, and generated plots.

## Details

This function is independent from
[`run_deggo()`](https://ymbouamboua.github.io/DEGgo/reference/run_deggo.md)
and operates on time-course expression data supplied as either a numeric
expression matrix or a DESeq2 `DESeqDataSet`.
