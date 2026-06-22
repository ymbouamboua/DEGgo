# Generate DEGgo PCA plot

Generate DEGgo PCA plot

## Usage

``` r
plot_pca(
  dds,
  metadata,
  output_dir = "DEGgo_out",
  filename = "PCA_Plot",
  intgroup = "condition",
  color_by = NULL,
  shape_by = NULL,
  title = "PCA Plot",
  width = 8,
  height = 6,
  dpi = 300,
  style = "classic"
)
```

## Arguments

- dds:

  DESeq2 dataset.

- metadata:

  Sample metadata.

- output_dir:

  Output directory.

- filename:

  Output filename without extension.

- intgroup:

  Metadata columns used for PCA coloring/grouping.

- color_by:

  Metadata column used for point color.

- shape_by:

  Metadata column used for point shape.

- title:

  Plot title.

- width:

  Plot width in inches.

- height:

  Plot height in inches.

- dpi:

  Plot resolution.

- style:

  Plot style passed to `.deggo_theme()`.

## Value

A list containing the VST object, PCA data, and plot.
