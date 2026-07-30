# Aggregate raw counts by an experimental unit

Aggregate raw counts by an experimental unit

## Usage

``` r
aggregate_counts_by_unit(
  counts,
  metadata,
  sample_col,
  unit_col,
  group_col,
  extra_cols = NULL
)
```

## Arguments

- counts:

  Gene-by-sample count matrix.

- metadata:

  Sample metadata.

- sample_col:

  Sample identifier column.

- unit_col:

  Experimental unit column, e.g. Cage.

- group_col:

  Group column.

- extra_cols:

  Optional columns required to be constant inside each unit.

## Value

List containing aggregated counts and metadata.
