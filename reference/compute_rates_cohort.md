# Wrapper function to compute statistics relating to incidence, prevalences and risks (ratios and differences) per group and between groups

Wrapper function to compute statistics relating to incidence,
prevalences and risks (ratios and differences) per group and between
groups

## Usage

``` r
compute_rates_cohort(
  aesifup_input,
  fupCol = "fup",
  pyrCol = "pyr",
  pyr_offset = "pyr_offset",
  eventCol = "eventCount",
  iptw = "ip_weight",
  target_aesi,
  type = "incidence",
  CImethod = "wilson",
  risk_type = "survival",
  CumulativeInc = TRUE,
  comparison_measures = TRUE,
  scale_IR = 10000,
  model_based_control = FALSE,
  weighted_IR = FALSE,
  end_risk,
  output_format = "data.table",
  minimum_count_for_comparative = 3
)
```

## Arguments

- aesifup_input:

  data.table, typically representing two groups (exposed and control),
  with columns specifying person-time, occurence of events, and group
  membership. Must contain a column `group` with values EXPOSED and
  CONTROL. May be outputted by `process_cohort_data`

- fupCol:

  name of column representing days of follow-up time. defaults to "fup"

- pyrCol:

  name of column representing person-years of follow-up time. defaults
  to "pyr"

- pyr_offset:

  name of column which is used as input to poisson (GEE) models if
  necessary, typically person-time. "pyr_offset"

- eventCol:

  name column identifying binary event count for this unit, "eventCount"

- iptw:

  name column containing (ipt or other) weights to be used in the
  analysis

- target_aesi:

  character string giving name of the outcome variable, for informative
  logging

- type:

  "incidence" or "prevalence", defining if incidence rates or prevalenes
  are to be computed

- CImethod:

  defaults to "wilson"; option passed to `est_inc_prev` defining formula
  for analytic CI computation

- risk_type:

  defaults "survival". if set to survival, function outputs 1-KM risk
  differences and ratios. Otherwise prevalence/incidence ratios and
  differences, as specified by `type`

- CumulativeInc:

  TRUE/FALSE, compute cumulative incidences using 1-KM. If not, skip and
  return -99 flag

- comparison_measures:

  controls whether 1-KM is based on the latest available follow-up in
  both exposed and control (TRUE) or just the maximum follow-up in the
  dataset as a whole (FALSE)

- scale_IR:

  what scaling to perform for incidence/prevalence, i.e., events per X.
  Defaults to 10000

- model_based_control:

  TRUE/FALSE. Use GEE to estimate incidence/prevalence in the control
  group (TRUE) or an analytic formula that assumes independent
  observations (FALSE)

- weighted_IR:

  use iptw column to weight the incidence/prevalence estimates per group
  (TRUE) or output unweighted incidence/prevalence in the exposed and
  control (FALSE)

- end_risk:

  end of the risk period (follow-up time at which cumulative incidence
  is evaluated); may be `NULL`

- output_format:

  defaults to `data.table`, otherwise returns data.frame.

- minimum_count_for_comparative:

  defaults to 3; should a minimum event count be applied in order to
  display results, all estimates relating to event counts less than this
  will be suppressed (flagged -77). Use `NULL` to disable

## Value

One-row `data.table` (or data.frame) with the outcome name, patient
counts, outcome counts, person-years and incidence/prevalence per group,
cumulative incidences at `end_risk`, and crude and adjusted (weighted)
risk ratios and differences with 95% CIs.

## Details

Flag values replace estimates that could not be computed: -88 (model
failure or no data), -99 (not requested, e.g. `CumulativeInc = FALSE`),
-77 (suppressed by `minimum_count_for_comparative`).
