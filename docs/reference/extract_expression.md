# Extract normalized gene expression

Extract normalized or variance-stabilized expression values from a
DESeq2 object and return a long-format table.

## Usage

``` r
extract_expression(
  dds,
  metadata,
  genes,
  assay = c("vst", "normalized", "log2_normalized", "raw"),
  gene_col = "SYMBOL",
  collapse_fun = sum
)
```

## Arguments

- dds:

  DESeqDataSet object.

- metadata:

  Sample metadata.

- genes:

  Character vector of gene symbols or identifiers.

- assay:

  One of `"vst"`, `"normalized"`, `"log2_normalized"` or `"raw"`.

- gene_col:

  Character. Column containing gene names or identifiers in
  rowData(dds).

- collapse_fun:

  Function used to collapse duplicated gene rows.

## Value

Long-format expression table.
