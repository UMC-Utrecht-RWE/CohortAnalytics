# function to compute incidence/prevalence rate ratios from a GEE model fit object
# also computes incidence/prevalence rate differences + CIs with the delta method


# takes a model object
# extracts point estimate and standard error
# computes CI
# exponentiates to get incidence rate / prevalence proportion ratios, CIs
# extract IR/PP Difference and compute SEs as appropriate (msm::deltamethod)
# note the same computation happens on models for incidence or prevalence rates
# depends on msm ; if analytic calcualtion on, also fmsb

extract_ests_gee <- function(
  modelobj, # geeglm or bigglm output model object
  scale_IR = 100000, # how to scale IR differences
  throw_warning = FALSE
) {
  if (is.character(modelobj)) {
    if (throw_warning) {
      ("warning: model object is character string, returning NA flags -88")
    }
    return(data.frame(
      irr_est = -88,
      irr_lb = -88,
      irr_ub = -88,
      ird_est = -88,
      ird_lb = -88,
      ird_ub = -88
    ))
  }

  # Detect model type and extract parameters accordingly
  if (inherits(modelobj, "bigglm")) {
    # Extract from bigglm object: summary()$mat has columns Coef, SE, p
    model_summary <- summary(modelobj)
    mat <- model_summary$mat

    int_par <- mat[1, "Coef"]
    int_se <- mat[1, "SE"]
    group_par <- mat[2, "Coef"]
    group_se <- mat[2, "SE"]

    # Get variance-covariance matrix via vcov()
    var_matrix <- vcov(modelobj)
    coef_vector <- coef(modelobj)
  } else {
    # Assume geeM object: use beta and se.robust from summary
    model_summary <- summary(modelobj)
    int_par <- model_summary$beta[1]
    int_se <- model_summary$se.robust[1]
    group_par <- model_summary$beta[2]
    group_se <- model_summary$se.robust[2]

    var_matrix <- modelobj$var
    coef_vector <- coef(modelobj)
  }

  # collect parameters
  par_ests <- data.frame(
    irr_est = group_par,
    irr_lb = group_par + qnorm(0.025) * group_se,
    irr_ub = group_par + qnorm(0.975) * group_se
  )

  # take exponential to get rate ratios
  rr_out <- exp(par_ests)

  # ---- calculate risk differences using delta method ----
  rd_est <- exp(int_par + group_par) - exp(int_par)
  rd_se <- msm::deltamethod(~ exp(x1 + x2) - exp(x1), coef_vector, var_matrix)
  rd_out <- data.frame(
    rd_est = rd_est,
    rd_lb = rd_est + qnorm(0.025) * rd_se,
    rd_ub = rd_est + qnorm(0.975) * rd_se
  ) * scale_IR

  # save and process output object
  outobj <- data.frame(rr_out, rd_out)
  rownames(outobj) <- NULL
  return(outobj)
}
