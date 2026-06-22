# Generate a PowerPoint report from DEGgo results

Generate a PowerPoint report from DEGgo results

## Usage

``` r
generate_deggo_pptx(
  results,
  output_file = "DEGgo_Report.pptx",
  title = "DEGgo RNA-seq report",
  subtitle = "Differential expression, visualization and Gene Ontology summary"
)
```

## Arguments

- results:

  DEGgo result object returned by run_deggo().

- output_file:

  Output PPTX file.

- title:

  Report title.

- subtitle:

  Report subtitle.

## Value

Path to the generated PPTX file.
