# Incidence rate or prevalence with analytic confidence intervals

Incidence uses exact Poisson (chi-squared) CIs; prevalence uses Wilson
or Clopper-Pearson CIs. CIs assume independent observations.

## Usage

``` r
est_inc_prev(
  n_pat,
  n_out,
  py = NULL,
  scale_IR,
  type = "incidence",
  CImethod = NULL
)
```

## Arguments

- n_pat:

  number of patients (individuals)

- n_out:

  number of outcomes

- py:

  person years of follow up

- scale_IR:

  how should the incidence/prevalence be scaled; i.e., incidence per X
  units of person-time

- type:

  "incidence" or "prevalence"

- CImethod:

  "wilson" or "clopper" based estimation of CIs

## Value

data.frame with `ir_est`, `ir_lb`, `ir_ub`, scaled by `scale_IR`. All
-88 if `py == 0`.
