# Plot GO enrichment terms across all DEGgo comparisons

Generate Gene Ontology enrichment plots for all comparisons stored in a
DEGgo results object. This function iterates over `results$go_results`,
creates one GO plot per comparison using
[`plot_go_terms()`](https://ymbouamboua.github.io/DEGgo/reference/plot_go_terms.md),
and optionally saves each plot as PNG and PDF.

## Usage

``` r
plot_all_go_terms(
  results,
  top_n = 10,
  txtsize = 8,
  style = "bw",
  output_dir = NULL,
  width = 8,
  height = 6,
  dpi = 300
)
```

## Arguments

- results:

  A DEGgo results object returned by
  [`run_deggo()`](https://ymbouamboua.github.io/DEGgo/reference/run_deggo.md).
  Must contain a `go_results` list.

- top_n:

  Integer. Number of top enriched GO terms to display per plot. Default
  is `10`.

- txtsize:

  Numeric. Base text size used in the GO plots. Default is `8`.

- style:

  Character. Plot theme style passed to
  [`plot_go_terms()`](https://ymbouamboua.github.io/DEGgo/reference/plot_go_terms.md).
  Default is `"bw"`.

- output_dir:

  Optional character. Directory where GO plots are saved. If `NULL`,
  plots are returned but not written to disk.

- width:

  Numeric. Plot width in inches for saved files. Default is `8`.

- height:

  Numeric. Plot height in inches for saved files. Default is `6`.

- dpi:

  Integer. Resolution in dots per inch for saved PNG files. Default is
  `300`.

## Value

A named list of `ggplot` objects, one per comparison. Comparisons
without enriched GO terms return `NULL`.

## Examples

``` r
if (FALSE) { # \dontrun{
go_plots <- plot_all_go_terms(
  results = results,
  top_n = 10,
  output_dir = "GO_plots"
)
} # }
```
