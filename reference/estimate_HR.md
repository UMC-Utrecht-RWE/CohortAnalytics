# Estimate a hazard ratio for exposed vs control

Fits a Cox model (`group` as the only covariate) with robust standard
errors clustered on `person_id`.

## Usage

``` r
estimate_HR(
  aesifup_input,
  fupCol = "fup",
  eventCol = "eventCount",
  iptw = "ip_weights",
  model_type = "crude",
  return_dummy_output = FALSE,
  aesi_name = ""
)
```

## Arguments

- aesifup_input:

  data.table with `group`, `person_id` and the follow-up, event and
  weight columns.

- fupCol:

  name of the follow-up time column.

- eventCol:

  name of the event indicator column.

- iptw:

  name of the weight column, used when `model_type = "adj"`.

- model_type:

  "crude" (unweighted) or "adj" (weighted by `iptw`).

- return_dummy_output:

  if TRUE, return the flag output without fitting.

- aesi_name:

  outcome name, for logging.

## Value

data.frame with `hr_est`, `hr_lb`, `hr_ub` (95% CI). All -88 if the
model cannot be fitted or there are no events.
