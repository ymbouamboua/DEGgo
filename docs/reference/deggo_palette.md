# Retrieve a DEGgo color palette

Returns a qualitative, sequential, or diverging color palette suitable
for DEGgo visualizations and annotation variables.

## Usage

``` r
deggo_palette(
  palette = c("default", "nature", "jama", "nejm", "lancet", "viridis", "okabe"),
  n = NULL,
  type = c("discrete", "sequential", "diverging"),
  direction = 1
)
```

## Arguments

- palette:

  Character string specifying the palette name. Available palettes are
  `"default"`, `"nature"`, `"jama"`, `"nejm"`, `"lancet"`, `"viridis"`,
  and `"okabe"`.

- n:

  Optional positive integer giving the number of colors to return. When
  `NULL`, the complete palette is returned.

- type:

  Character string specifying the palette type. One of `"discrete"`,
  `"sequential"`, or `"diverging"`.

- direction:

  Integer equal to `1` or `-1`. Use `-1` to reverse the palette.

## Value

A character vector of hexadecimal colors.

## Examples

``` r
deggo_palette("default")
#>  [1] "#0072B2" "#D55E00" "#009E73" "#CC79A7" "#E69F00" "#56B4E9" "#882255"
#>  [8] "#44AA99" "#AA4499" "#999999" "#000000" "#F0E442"
deggo_palette("okabe", n = 4)
#> [1] "#E69F00" "#56B4E9" "#009E73" "#F0E442"
deggo_palette("viridis", n = 8, type = "sequential")
#> [1] "#440154" "#482878" "#3E4989" "#31688E" "#26828E" "#1F9E89" "#35B779"
#> [8] "#6DCD59"
deggo_palette("default", n = 11, type = "diverging")
#>  [1] "#2166AC" "#3C8ABE" "#72B1D3" "#ABD1E5" "#D8E8F1" "#F7F7F7" "#FBE0D0"
#>  [8] "#F7BA9D" "#E8896C" "#CE5146" "#B2182B"
```
