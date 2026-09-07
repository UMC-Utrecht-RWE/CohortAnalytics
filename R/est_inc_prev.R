
# helper functions to compute incidence rates and prevalence rates per group

# this function uses analytic/closed form expressions (CIs valid under assumption of independent observations)
#' Estimate incidence rates or prevalence rates
#'
#' @param n_pat number of patients (individuals)
#' @param n_out number of outcomes
#' @param py person years of follow up
#' @param scale_IR how should the incidence/prevalence be scaled; i.e., incidence per X units of person-time
#' @param type "incidence" or "prevalence"
#' @param CImethod "wilson" or "clopper" or "bootstrap" or NULL
#' @param boot_vec vector of bootstrap result to pass to CImethod = "bootstrap"
#'
#' @export
est_inc_prev <- function(
    n_pat,
    n_out,
    py = NULL,
    scale_IR,
    type = "incidence",
    CImethod = NULL,
    boot_vec = NULL,
    py_zero_code = -88){
  #CImethod required for prevalence CIs
  # n_pat = n_pat_exp
  # n_out = n_out_exp
  # py = py_exp
  # scale_IR = scale_IR
  # type = type
  if(py == 0){ return(data.frame(ir = py_zero_code,lb_ir = py_zero_code, ub_ir = py_zero_code)) }
  if(type == "incidence"){
    if(is.null(py)){stop("for incidence calculation, py must be non-NULL")}
    # incidence CI using dobson formula
    ir <- n_out/py
    lb_ir <- max(0,((qchisq(0.025, df=(2*n_out)))/2)/py)
    ub_ir <- (qchisq(0.975, df=(2*(n_out+1)))/2)/py
  }

  if(type == "prevalence"){
    if (is.null(CImethod)) {
      ir <- n_out / n_pat
      lb_ir <- ub_ir <- NA
    } else if(CImethod == "clopper") {
      # prevalence using clopper formula
      test_out <- binom.test(n_out, n_pat)
      lb_ir <- test_out$conf.int[1]
      ub_ir <- test_out$conf.int[2]
      ir <- test_out$estimate
    } else if(CImethod == "wilson"){
      # prevalence using wilson formula
      ir <- n_out / n_pat
      z <- qnorm(0.975)  # z-score for 95% CI
      denom <- 1 + (z^2 / n_pat)
      term1 <- ir + (z^2 / (2 * n_pat))
      term2 <- z * sqrt((ir * (1 - ir) / n_pat) + (z^2 / (4 * n_pat^2)))

      lb_wilson <- (term1 - term2) / denom
      ub_wilson <- (term1 + term2) / denom

      # Use Wilson CI as the main output
      lb_ir <- lb_wilson
      ub_ir <- ub_wilson
    } else if (CImethod == "bootstrap") {
      # CI using 2.5 and 97.5 percentiles of bootstrap result
      # Validate boot_vec
      if (is.null(boot_vec) || length(boot_vec) == 0) {
        stop("`boot_vec` must be a non-empty numeric vector for bootstrap CIs.")
      }
      if (!is.numeric(boot_vec)) {
        stop("`boot_vec` must be numeric.")
      }

      # Compute percentile CI (2.5% and 97.5%)
      qs <- stats::quantile(boot_vec, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)

      lb_ir <- qs[1]
      ub_ir <- qs[2]
      ir <- NA
    }
  }

  return(data.table::data.table(ir_est = ir*scale_IR, ir_lb = lb_ir*scale_IR, ir_ub = ub_ir*scale_IR))
}


# model based estimation of ir's
# takes as input a model object (typicaly a GEE model)
est_inc_prev_model <- function(modelobj, group, scale_IR){
  if(is.character(modelobj)){ return(data.frame(ir_est = -88,ir_lb = -88, ir_ub = -88)) }

  coef_vec <- tryCatch(stats::coef(modelobj), error = function(e) NULL)
  if (is.null(coef_vec)) {
    return(data.frame(ir_est = -88, ir_lb = -88, ir_ub = -88))
  }

  int_idx <- which(names(coef_vec) == "(Intercept)")
  if (length(int_idx) == 0) {
    int_idx <- 1
  }
  int_idx <- int_idx[1]

  vcov_mat <- tryCatch(stats::vcov(modelobj), error = function(e) NULL)
  if (is.null(vcov_mat) && !is.null(modelobj$var)) {
    vcov_mat <- modelobj$var
  }
  if (is.null(vcov_mat)) {
    return(data.frame(ir_est = -88, ir_lb = -88, ir_ub = -88))
  }

  int_par <- coef_vec[int_idx]
  int_se <- sqrt(vcov_mat[int_idx, int_idx])

  # collect parameters
  # For one-group model fits, the intercept is the group-specific estimate.
  ir <- exp(int_par)
  lb_ir <- exp(int_par + qnorm(0.025)*int_se)
  ub_ir <- exp(int_par + qnorm(0.975)*int_se)

  return(data.frame(ir_est = ir*scale_IR,
                    ir_lb = lb_ir*scale_IR,
                    ir_ub =  ub_ir*scale_IR
  )
  )}#end of function
