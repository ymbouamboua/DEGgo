# Plot GO enrichment terms for all pairwise comparisons

Generates GO enrichment plots for every comparison returned by
[`run_deggo()`](https://ymbouamboua.github.io/DEGgo/reference/run_deggo.md)
and optionally exports figures to disk.

## Usage

``` r
plot_all_go_terms(
  results,
  top_n = 10,
  font_size = 8,
  output_dir = NULL,
  width = 8,
  height = 6,
  dpi = 300
)
```

## Arguments

- results:

  DEGgo results object returned by
  [`run_deggo()`](https://ymbouamboua.github.io/DEGgo/reference/run_deggo.md).

- top_n:

  Number of GO terms to display per regulation class (Up and Down).
  Default is `10`.

- font_size:

  Base font size used in plots. Default is `8`.

- output_dir:

  Optional output directory for exported figures. If `NULL`, plots are
  returned but not saved.

- width:

  Plot width in inches. Default is `8`.

- height:

  Plot height in inches. Default is `6`.

- dpi:

  Plot resolution. Default is `300`.

## Value

A named list of ggplot objects, one per comparison.

## Details

This function iterates over all pairwise comparisons stored in
`results$go_results` and calls
[`plot_go_terms()`](https://ymbouamboua.github.io/DEGgo/reference/plot_go_terms.md)
on each GO enrichment table.

If `output_dir` is supplied, PNG and PDF versions of all plots are
automatically exported.

## Examples

``` r
if (FALSE) { # \dontrun{

go_plots <- plot_all_go_terms(
  results = results,
  top_n = 10,
  font_size = 8,
  output_dir = "GO_plots"
)

names(go_plots)

go_plots$WAT_Female_PAMH_vs_PBS

} # }
```
