# Generate a DEGgo report

Generate a DEGgo report

## Usage

``` r
generate_deggo_report(
  results,
  output_dir,
  formats = "html",
  report_template = NULL,
  project_name = NULL
)
```

## Arguments

- results:

  DEGgo results object.

- output_dir:

  Output directory.

- formats:

  Report formats: `"html"`, `"pdf"`, or both.

- report_template:

  Optional custom R Markdown template.

- project_name:

  Optional project name.

## Value

Character vector containing generated report paths.
