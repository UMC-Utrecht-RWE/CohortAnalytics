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
  py_zero_code = -88
) {
  # CImethod required for prevalence CIs
  # n_pat = n_pat_exp
  # n_out = n_out_exp
  # py = py_exp
  # scale_IR = scale_IR
  # type = type
  if (py == 0) {
    return(data.frame(ir = py_zero_code, lb_ir = py_zero_code, ub_ir = py_zero_code))
  }
  if (type == "incidence") {
    if (is.null(py)) {
      stop("for incidence calculation, py must be non-NULL")
    }
    # incidence CI using dobson formula
    ir <- n_out / py
    lb_ir <- max(0, ((qchisq(0.025, df = (2 * n_out))) / 2) / py)
    ub_ir <- (qchisq(0.975, df = (2 * (n_out + 1))) / 2) / py
  }

  if (type == "prevalence") {
    if (is.null(CImethod)) {
      ir <- n_out / n_pat
      lb_ir <- ub_ir <- NA
    } else if (CImethod == "clopper") {
      # prevalence using clopper formula
      test_out <- binom.test(n_out, n_pat)
      lb_ir <- test_out$conf.int[1]
      ub_ir <- test_out$conf.int[2]
      ir <- test_out$estimate
    } else if (CImethod == "wilson") {
      # prevalence using wilson formula
      ir <- n_out / n_pat
      z <- qnorm(0.975) # z-score for 95% CI
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

  return(data.table::data.table(ir_est = ir * scale_IR, ir_lb = lb_ir * scale_IR, ir_ub = ub_ir * scale_IR))
}


# model based estimation of ir's
# takes as input a model object (typicaly a GEE model)
#' Estimate model-based incidence/prevalence for a specific group from a GEE Poisson model
#'
#' Extracts a model-based incidence/prevalence estimate and Wald-type 95% CI on the
#' log scale, then exponentiates and scales the result.
#'
#' Supports two model structures:
#' 1) subgroup-specific model without a group term (intercept-only for that subgroup), and
#' 2) model with a group term, where CONTROL is based on the intercept and EXPOSED is
#' based on intercept + group coefficient.
#'
#' For EXPOSED in a model with a group term, the variance of the linear predictor is
#' computed as var(intercept) + var(group) + 2*cov(intercept, group).
#'
#' @param modelobj fitted model object (typically `geeM::geem`) or a character error flag
#' @param group target group label; must be `"CONTROL"` or `"EXPOSED"`
#' @param scale_IR multiplicative scaling factor for the output estimate and CI bounds
#'
#' @return data.frame with `ir_est`, `ir_lb`, and `ir_ub`
est_inc_prev_model <- function(modelobj, group, scale_IR) {
  # Upstream code may return model-fit failures as character strings.
  # Preserve existing sentinel behavior for downstream pipelines.
  if (is.character(modelobj)) {
    return(data.frame(ir_est = -88, ir_lb = -88, ir_ub = -88))
  }

  # Restrict to the two supported group labels to avoid silent mis-specification.
  if (!(group %in% c("CONTROL", "EXPOSED"))) {
    stop("group must be either 'CONTROL' or 'EXPOSED'")
  }

  # Extract coefficients from the fitted model.
  model_summary <- summary(modelobj)
  beta <- stats::coef(modelobj)
  if (is.null(names(beta))) {
    names(beta) <- paste0("beta_", seq_along(beta))
  }

  # Try to use robust covariance when available; otherwise approximate using robust SEs.
  vcov_robust <- NULL
  if (!is.null(model_summary$cov.robust) && is.matrix(model_summary$cov.robust)) {
    vcov_robust <- model_summary$cov.robust
  } else if (!is.null(model_summary$vbeta.robust) && is.matrix(model_summary$vbeta.robust)) {
    vcov_robust <- model_summary$vbeta.robust
  } else if (!is.null(modelobj$vbeta) && is.matrix(modelobj$vbeta)) {
    vcov_robust <- modelobj$vbeta
  } else {
    vcov_robust <- tryCatch(stats::vcov(modelobj), error = function(e) NULL)
  }

  # Locate intercept and group-term positions in the coefficient vector.
  int_idx <- which(names(beta) == "(Intercept)")
  if (length(int_idx) == 0) {
    int_idx <- 1
  }

  group_idx <- grep("^group", names(beta))

  # If the model was fit within one subgroup only (no group term),
  # the intercept corresponds directly to that subgroup's incidence/prevalence.
  if (length(group_idx) == 0) {
    # Linear predictor for the subgroup-specific fit.
    eta <- unname(beta[int_idx])
    if (!is.null(vcov_robust) && nrow(vcov_robust) >= int_idx) {
      se_eta <- sqrt(vcov_robust[int_idx, int_idx])
    } else if (!is.null(model_summary$se.robust) && length(model_summary$se.robust) >= int_idx) {
      se_eta <- model_summary$se.robust[int_idx]
    } else {
      stop("Unable to extract standard error from model object")
    }
  } else {
    # Model includes a group indicator. CONTROL uses intercept only;
    # EXPOSED uses intercept + group coefficient.
    if (group == "CONTROL") {
      # Linear predictor for reference group.
      eta <- unname(beta[int_idx])
      if (!is.null(vcov_robust) && nrow(vcov_robust) >= int_idx) {
        se_eta <- sqrt(vcov_robust[int_idx, int_idx])
      } else if (!is.null(model_summary$se.robust) && length(model_summary$se.robust) >= int_idx) {
        se_eta <- model_summary$se.robust[int_idx]
      } else {
        stop("Unable to extract standard error from model object")
      }
    } else {
      grp_idx <- group_idx[1]
      # Linear predictor for EXPOSED: intercept + group effect.
      eta <- unname(beta[int_idx] + beta[grp_idx])

      if (!is.null(vcov_robust) &&
        nrow(vcov_robust) >= max(int_idx, grp_idx) &&
        ncol(vcov_robust) >= max(int_idx, grp_idx)) {
        # Delta-method variance for a sum of coefficients.
        var_eta <- vcov_robust[int_idx, int_idx] +
          vcov_robust[grp_idx, grp_idx] +
          2 * vcov_robust[int_idx, grp_idx]
        se_eta <- sqrt(max(0, var_eta))
      } else if (!is.null(model_summary$se.robust) && length(model_summary$se.robust) >= max(int_idx, grp_idx)) {
        # Fallback approximation if covariance is unavailable.
        se_eta <- sqrt(model_summary$se.robust[int_idx]^2 + model_summary$se.robust[grp_idx]^2)
      } else {
        stop("Unable to extract standard error from model object")
      }
    }
  }

  # Convert back from log scale and form 95% Wald CI.
  ir <- exp(eta)
  lb_ir <- exp(eta + qnorm(0.025) * se_eta)
  ub_ir <- exp(eta + qnorm(0.975) * se_eta)

  return(data.frame(
    ir_est = ir * scale_IR,
    ir_lb = lb_ir * scale_IR,
    ir_ub = ub_ir * scale_IR
  ))
} # end of function
