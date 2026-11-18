# Pre-process cohort data for analysis

Standard pre-processing before
[`compute_rates_cohort()`](https://umc-utrecht-rwe.github.io/CohortAnalytics/reference/compute_rates_cohort.md):
sets `CONTROL` as the reference level of `group`, creates a numeric
person id (`person_id_num`), sorts by it, creates the log person-time
offset (`pyr_offset`), and drops pairs with zero follow-up or missing
weights.

## Usage

``` r
process_cohort_data(
  aesifup,
  estimator,
  fupdays = "fup",
  pyrCol = "pyr",
  ipweight = "ip_weight"
)
```

## Arguments

- aesifup:

  data.table with columns `group` (EXPOSED/CONTROL), `person_id`, `id`
  (pair id) and follow-up. Modified by reference for some steps.

- estimator:

  "incidence" or "prevalence". For prevalence the offset is `log(1)`.

- fupdays:

  name of the follow-up (days) column.

- pyrCol:

  name of the person-years column; computed as `fup / 365.25` if absent.

- ipweight:

  name of the weight column; rows with `NA` weights are dropped with
  their pair.

## Value

data.table with added columns `person_id_num` and `pyr_offset`.
