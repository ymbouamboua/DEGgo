# Extract genes found in GO terms

Extract genes found in GO terms

## Usage

``` r
deggo_extract_go_genes(go_results, genes, gene_sep = "/")
```

## Arguments

- go_results:

  GO enrichment result data.frame.

- genes:

  Character vector of genes to search.

- gene_sep:

  Separator used in the geneID column.

## Value

A data.frame with GO terms containing selected genes.
