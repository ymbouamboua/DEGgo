# Prepare one GSE98965 tissue for circadian analysis

Extract one tissue from the GSE98965 baboon transcriptomic atlas,
construct a numeric genes-by-samples expression matrix, filter lowly
expressed genes, generate numeric Zeitgeber Time metadata, and prepare
gene-symbol annotations.

## Usage

``` r
prepare_circadian_tissue(
  expression_data,
  tissue_code = "LIV",
  min_expression = 1,
  min_samples = 3,
  log_transform = TRUE
)
```

## Arguments

- expression_data:

  A `data.frame` containing the processed GSE98965 FPKM expression
  table, usually returned by
  [`download_gse98965()`](https://ymbouamboua.github.io/DEGgo/reference/download_gse98965.md).

- tissue_code:

  Character string specifying the GSE98965 tissue code to extract. The
  default is `"LIV"` for liver.

- min_expression:

  Numeric value. Minimum FPKM expression required for a gene to be
  considered expressed in a sample.

- min_samples:

  Positive integer. Minimum number of samples in which a gene must reach
  `min_expression`.

- log_transform:

  Logical. If `TRUE`, expression values are transformed as
  `log2(FPKM + 1)`.

## Value

A named list containing:

- expr:

  Numeric genes-by-samples expression matrix.

- metadata:

  Sample metadata with `sample`, `tissue`, `ZT`, and numeric `time`
  columns.

- gene_annotation:

  Gene identifier and symbol annotation table.

- filter_summary:

  Summary of expression filtering.

- tissue_code:

  Selected tissue code.

- assay:

  Expression representation returned by the function.

## Details

Sample columns are detected using the pattern `TISSUE.ZTXX`, for example
`LIV.ZT02`.

Duplicate gene identifiers are preserved using
[`base::make.unique()`](https://rdrr.io/r/base/make.unique.html). The
resulting unique identifiers are stored in the `matrix_gene_id`
annotation column and correspond exactly to the expression matrix row
names.

## Examples

``` r
if (FALSE) { # \dontrun{
baboon_expression <- download_gse98965()

liver <- prepare_circadian_tissue(
  expression_data = baboon_expression,
  tissue_code = "LIV",
  min_expression = 1,
  min_samples = 3
)

dim(liver$expr)
head(liver$metadata)
head(liver$gene_annotation)
} # }
```
