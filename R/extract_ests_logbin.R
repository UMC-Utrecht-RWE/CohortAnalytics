#' Extract risk ratio estimates from a log-binomial model
#'
#' Retrieves the point estimate and Wald 95% confidence interval for the
#' EXPOSED group coefficient from a `glm` object returned by
#' `fitmod_logbin`.  Returns dummy codes when model fitting failed.
#'
#' @param modelobj `glm` object from `fitmod_logbin`, or a character string
#'   if the model failed to converge
#' @param dummy_code scalar returned for each estimate when model is a failure
#'
#' @return a one-row `data.frame` with columns `rr_est`, `rr_lb`, `rr_ub`
#' @export
extract_ests_logbin <- function(modelobj, dummy_code = -88) {
  dummy_out <- data.frame(rr_est = dummy_code,
                          rr_lb  = dummy_code,
                          rr_ub  = dummy_code)

  if (is.character(modelobj)) return(dummy_out)

  tryCatch({
    coefs   <- coef(modelobj)
    vcovmat <- vcov(modelobj)

    # locate coefficient for EXPOSED (robust to any reference ordering)
    idx <- grep("EXPOSED", names(coefs))
    if (length(idx) == 0) return(dummy_out)

    log_rr <- coefs[idx]
    log_se <- sqrt(vcovmat[idx, idx])

    data.frame(
      rr_est = exp(log_rr),
      rr_lb  = exp(log_rr - qnorm(0.975) * log_se),
      rr_ub  = exp(log_rr + qnorm(0.975) * log_se)
    )
  }, error = function(e) dummy_out)
}
