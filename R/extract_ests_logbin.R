#' Extract risk ratio estimates from a log-binomial model
#'
#' Retrieves the point estimate and Wald 95% confidence interval for the
#' EXPOSED group coefficient from a `glm` object returned by
#' `fitmod_logbin`.  Returns dummy codes when model fitting failed.
#'
#' @param modelobj `glm` object from `fitmod_logbin`, or a character string
#'   if the model failed to converge
#' @param dummy_code scalar returned for each estimate when model is a failure
#' @param return_format either "rr" (default) to return `rr_est`, `rr_lb`,
#'   `rr_ub`, or "risk" to return `irr_*` and `ird_*` columns compatible with
#'   `extract_ests_gee` output.
#'
#' @return a one-row `data.frame`; column names depend on `return_format`.
#' @export
extract_ests_logbin <- function(modelobj,
                                dummy_code = -88,
                                return_format = c("rr", "risk")) {
  return_format <- match.arg(return_format)

  dummy_rr <- data.frame(rr_est = dummy_code,
                         rr_lb  = dummy_code,
                         rr_ub  = dummy_code)

  dummy_risk <- data.frame(irr_est = dummy_code,
                           irr_lb  = dummy_code,
                           irr_ub  = dummy_code,
                           ird_est = dummy_code,
                           ird_lb  = dummy_code,
                           ird_ub  = dummy_code)

  if (return_format == "risk") {
    base_out <- dummy_risk
  } else {
    base_out <- dummy_rr
  }

  if (is.character(modelobj)) return(base_out)

  tryCatch({
    coefs   <- coef(modelobj)
    vcovmat <- vcov(modelobj)

    # locate coefficient for EXPOSED (robust to any reference ordering)
    idx <- grep("EXPOSED", names(coefs))
    if (length(idx) == 0) return(base_out)

    log_rr <- coefs[idx]
    log_se <- sqrt(vcovmat[idx, idx])

    rr_out <- data.frame(
      rr_est = exp(log_rr),
      rr_lb  = exp(log_rr - qnorm(0.975) * log_se),
      rr_ub  = exp(log_rr + qnorm(0.975) * log_se)
    )

    if (return_format == "rr") {
      return(rr_out)
    }

    data.frame(irr_est = rr_out$rr_est,
               irr_lb  = rr_out$rr_lb,
               irr_ub  = rr_out$rr_ub,
               ird_est = -88,
               ird_lb  = -88,
               ird_ub  = -88)
  }, error = function(e) base_out)
}
