# Plot normalized gene expression

Generates publication-ready violin, boxplot, or barplot visualizations
from normalized expression values produced by
`extract_normalized_expression()`.

## Usage

``` r
plot_gene_expression(
  expr_df,
  gene,
  x = "treatment",
  color = x,
  facet = NULL,
  geom = c("violin", "boxplot", "barplot"),
  comparisons = list(c("PBS", "PAMH")),
  stats = TRUE,
  stat_method = "wilcox.test",
  stat_label = c("p.signif", "p.format"),
  output_dir = NULL,
  filename = NULL,
  width = 6,
  height = 5,
  dpi = 300,
  style = "classic",
  txtsize = 12,
  x_ang = 45,
  ncol = NULL
)
```

## Arguments

- expr_df:

  Long-format expression table.

- gene:

  Gene symbol to plot.

- x:

  Metadata variable displayed on the x-axis.

- color:

  Metadata variable used for fill colors.

- facet:

  Optional metadata variable used for faceting.

- geom:

  Plot type: `"violin"`, `"boxplot"`, or `"barplot"`.

- comparisons:

  List of comparisons passed to
  [`ggpubr::stat_compare_means()`](https://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html).

- stats:

  Logical. Add statistical comparisons.

- stat_method:

  Statistical test method.

- stat_label:

  Label style: `"p.signif"` or `"p.format"`.

- output_dir:

  Output directory.

- filename:

  Output filename.

- width:

  Plot width.

- height:

  Plot height.

- dpi:

  Plot resolution.

- style:

  Theme style passed to `.deggo_theme()`.

- txtsize:

  Base font size.

- x_ang:

  Numeric. Angle of x-axis labels.

- ncol:

  Number of plot columns.

## Value

A ggplot object.

## Details

Statistical comparisons can be added using
[`ggpubr::stat_compare_means()`](https://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html).

## Examples

``` r
if (FALSE) { # \dontrun{
plot_gene_expression(
  expr_df,
  gene = "Adipoq",
  x = "treatment",
  color = "treatment",
  facet = "tissue",
  geom = "violin"
)
} # }
```
