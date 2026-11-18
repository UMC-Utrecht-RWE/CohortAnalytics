#' Estimate a hazard ratio for exposed vs control
#'
#' Fits a Cox model (`group` as the only covariate) with robust standard errors
#' clustered on `person_id`.
#'
#' @param aesifup_input data.table with `group`, `person_id` and the follow-up, event and weight columns.
#' @param fupCol name of the follow-up time column.
#' @param eventCol name of the event indicator column.
#' @param iptw name of the weight column, used when `model_type = "adj"`.
#' @param model_type "crude" (unweighted) or "adj" (weighted by `iptw`).
#' @param return_dummy_output if TRUE, return the flag output without fitting.
#' @param aesi_name outcome name, for logging.
#'
#' @returns data.frame with `hr_est`, `hr_lb`, `hr_ub` (95% CI). All -88 if the
#'   model cannot be fitted or there are no events.
#' @export
estimate_HR <- function(
  aesifup_input,
  fupCol = "fup",
  eventCol = "eventCount",
  iptw = "ip_weights",
  model_type = "crude",
  return_dummy_output = FALSE,
  aesi_name = ""
) {
  dummy_output <- data.frame(hr_est = -88, hr_lb = -88, hr_ub = -88)
  if (return_dummy_output == TRUE) {
    return(dummy_output)
  }

  # check if non-zero events in both groups, return dummy output
  eventsums <- aesifup_input[, sum(get(eventCol)), by = group]
  if (!any(eventsums[, 2] > 0)) {
    zero_groups <- paste0(
      as.character(eventsums[V1 == 0]$group),
      collapse = ","
    )
    logger::log_info(paste0(
      "HR model estimation fails, zero events in groups ",
      zero_groups
    ))
    return(dummy_output)
  }

  # fit coxmodel
  model_formula <- stats::as.formula(paste0(
    "survival::Surv(",
    fupCol,
    ",",
    eventCol,
    ") ~ group"
  ))
  modelobj <- fitmod_cox(
    model_formula,
    iptw = iptw,
    model_type = model_type,
    aesi_name = aesi_name,
    aesifup_input
  )

  # if model fitting error, return error flag

  if (is.character(modelobj)) {
    return(dummy_output)
  }

  # get coefficients
  coxcoef <- summary(modelobj)$coefficients
  hr_est <- coxcoef[colnames(coxcoef) == "coef"]
  hr_se <- coxcoef[colnames(coxcoef) == "robust se"]

  # get confidence interval
  hr_lb <- hr_est + stats::qnorm(0.025) * hr_se
  hr_ub <- hr_est + stats::qnorm(0.975) * hr_se

  # exponentiate
  hr_est <- exp(hr_est)
  hr_lb <- exp(hr_lb)
  hr_ub <- exp(hr_ub)

  # return object
  return(data.frame(hr_est = hr_est, hr_lb = hr_lb, hr_ub = hr_ub))
}
