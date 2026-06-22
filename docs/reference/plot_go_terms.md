# Plot GO terms by regulation status

Plot GO terms by regulation status

## Usage

``` r
plot_go_terms(
  go_df,
  go_terms_of_interest = NULL,
  comparison = NULL,
  top_n = 10,
  wrap_width = 50,
  color_values = c(Up = "#740001", Down = "#6497b1"),
  size_range = c(2, 6),
  style = "bw",
  x_angle = 0,
  txtsize = 12
)
```

## Arguments

- go_df:

  GO enrichment data frame from
  [`run_go_enrichment()`](https://ymbouamboua.github.io/DEGgo/reference/run_go_enrichment.md).

- go_terms_of_interest:

  Optional GO terms to display.

- comparison:

  Optional comparison name used as plot title.

- top_n:

  Number of GO terms to show per regulation.

- wrap_width:

  Width for wrapping GO term labels.

- color_values:

  Named colors for Up and Down.

- size_range:

  Point size range.

- style:

  Theme style passed to `.deggo_theme()`.

- x_angle:

  X-axis text angle.

- txtsize:

  Base font size.

## Value

A ggplot object.
