# ---- function to estimate hazard ratios -------------
# returns dummy output if not possible to estimate
#' Estimate hazard ratios using Cox proportional hazards model
#'
#' @param aesifup_input data
#' @param fupCol column name of follow up time, to be used in survival model
#' @param eventCol event column name
#' @param groupCol column identifying exposed/control group
#' @param iptw column name of inverse probability weight
#' @param model_type weighted or crude
#' @param return_dummy_output TRUE/FALSE whether to return a dummy output
#' @param aesi_name outcome label
#' @param dummy_code use code instead of NA
#' @export

estimate_HR <- function(aesifup_input, fupCol = "fup",
                        eventCol = "eventCount",
                        groupCol = "group",
                        iptw = "ip_weights",
                        model_type = "crude",
                        return_dummy_output = FALSE,
                        aesi_name = "",
                        dummy_code = -88) {
  dummy_output <- data.frame(hr_est = dummy_code, hr_lb = dummy_code, hr_ub = dummy_code)
  if (return_dummy_output == TRUE) {
    return(dummy_output)
  }

  # check if non-zero events in both groups, return dummy output

dt <- data.table::as.data.table(aesifup_input)

# Fully explicit column access (no get(), no .SD, no .())
eventsums <- data.table::data.table(
  group_val = dt[[groupCol]],
  event_val = dt[[eventCol]]
)[
  ,
  .(event_sum = sum(event_val, na.rm = TRUE)),
  by = group_val
]

data.table::setnames(eventsums, "group_val", groupCol)

if (nrow(eventsums) == 0 || !any(eventsums$event_sum > 0)) {
    zero_groups <- paste0(as.character(eventsums[event_sum == 0][[groupCol]]), collapse = ",")
    logger::log_info(paste0("HR model estimation fails, zero events in groups ", zero_groups))
    return(dummy_output)
  }

  # fit coxmodel
  model_formula <- as.formula(paste0("survival::Surv(", fupCol, ",", eventCol, ") ~ ", groupCol))
  modelobj <- fitmod_cox(model_formula,
    iptw = iptw,
    model_type = model_type, aesi_name = aesi_name, aesifup_input
  )


  # if model fitting returned logger object or not a coxph, return dummy output
  if (!inherits(modelobj, "coxph")) {
    return(dummy_output)
  }


  # if model fitting error, return error flag

  if (is.character(modelobj)) {
    return(dummy_output)
  }

  # get coefficients
  coxcoef <- summary(modelobj)$coefficients
  hr_est <- coxcoef[colnames(coxcoef) == "coef"]
  hr_se <- coxcoef[colnames(coxcoef) == "robust se"]
  if (length(hr_se) == 0) {
    hr_se <- 0
  }

  # get confidence interval
  hr_lb <- hr_est + qnorm(0.025) * hr_se
  hr_ub <- hr_est + qnorm(0.975) * hr_se

  # exponentiate
  hr_est <- exp(hr_est)
  hr_lb <- exp(hr_lb)
  hr_ub <- exp(hr_ub)

  # return object
  return(data.frame(hr_est = hr_est, hr_lb = hr_lb, hr_ub = hr_ub))
}
