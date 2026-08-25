#' Summarize predicted probabilities
#'
#' Computes the sample size, mean, standard deviation, and Wald-style
#' confidence interval for a vector of predicted probabilities.
#'
#' @param pred_values Numeric vector of predicted probabilities.
#' @param conf_level Confidence level for the Wald interval.
#'
#' @return A named list with `n`, `mean`, `sd`, `ci_lb`, and `ci_ub`.
#' @keywords internal
.summarize_timevarying_predictions <- function(pred_values,
                                               conf_level = 0.95) {
  pred_values <- pred_values[!is.na(pred_values)]
  n_pred <- length(pred_values)

  if (n_pred == 0) {
    return(list(n = 0L,
                mean = NA_real_,
                sd = NA_real_,
                ci_lb = NA_real_,
                ci_ub = NA_real_))
  }

  pred_mean <- mean(pred_values)
  pred_sd <- stats::sd(pred_values)

  if (n_pred <= 1 || is.na(pred_sd)) {
    return(list(n = n_pred,
                mean = pred_mean,
                sd = pred_sd,
                ci_lb = NA_real_,
                ci_ub = NA_real_))
  }

  z_value <- stats::qnorm(1 - (1 - conf_level) / 2)
  pred_se <- pred_sd / sqrt(n_pred)

  list(n = n_pred,
       mean = pred_mean,
       sd = pred_sd,
       ci_lb = pred_mean - z_value * pred_se,
       ci_ub = pred_mean + z_value * pred_se)
}