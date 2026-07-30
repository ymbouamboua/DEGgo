# Print a DEGgo design quality control summary

Displays a concise summary of the study design assessment returned by
[`design_qc()`](https://ymbouamboua.github.io/DEGgo/reference/design_qc.md).
The summary includes the number of samples and groups, whether repeated
measurements are present, the dominant experimental factor, the
recommended differential expression framework, and any detected design
warnings.

## Usage

``` r
# S3 method for class 'deggo_design_qc'
print(x, ...)
```

## Arguments

- x:

  An object of class `"deggo_design_qc"` returned by
  [`design_qc()`](https://ymbouamboua.github.io/DEGgo/reference/design_qc.md).

- ...:

  Additional arguments passed to or from other methods.

## Value

The input object, invisibly.

## See also

[`design_qc()`](https://ymbouamboua.github.io/DEGgo/reference/design_qc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- design_qc(
  counts = counts,
  metadata = metadata,
  sample_col = "sample",
  group_col = "group"
)

print(res)
} # }
```
