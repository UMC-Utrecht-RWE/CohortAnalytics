#' Predict monthly treated and untreated probabilities
#'
#' Builds a month-level prediction grid and obtains one predicted probability
#' under treated and untreated status for each observed time value.
#'
#' @param fit A fitted time-varying model object.
#' @param model_data Data used as the prediction template.
#' @param timeVar Name of the time index column.
#' @param idCol Name of the person identifier column.
#' @param group_col Name of the treatment group column.
#' @param treated_group Value representing the treated group.
#' @param untreated_group Value representing the untreated group.
#'
#' @return A data.table with one row per observed time value and summary
#'   statistics for the treated and untreated predicted probabilities.
#' @keywords internal
.predict_timevarying_probabilities <- function(
    fit,
    model_data,
    timeVar,
    idCol = "person_id_num",
    group_col = "group",
    treated_group = "EXPOSED",
    untreated_group = "CONTROL",
    conf_level = 0.95) {

  if (is.character(fit)) {
    return(.empty_timevarying_prediction_table(timeVar = timeVar))
  }

  df <- as.data.frame(model_data)

  if (is.null(timeVar) || !(timeVar %in% names(df))) {
    stop("A valid timeVar column is required for time-varying predictions")
  }

  if (!(group_col %in% names(df))) {
    stop(paste0("Column '", group_col, "' not found in model_data"))
  }

  month_values <- sort(unique(df[[timeVar]][!is.na(df[[timeVar]])]))

  if (length(month_values) == 0) {
    return(.empty_timevarying_prediction_table(timeVar = timeVar))
  }

  # Number of observations per month and treatment group
  treated_counts <- vapply(
    month_values,
    function(month_value) {
      sum(
        df[[timeVar]] == month_value &
          df[[group_col]] == treated_group,
        na.rm = TRUE
      )
    },
    integer(1)
  )

  untreated_counts <- vapply(
    month_values,
    function(month_value) {
      sum(
        df[[timeVar]] == month_value &
          df[[group_col]] == untreated_group,
        na.rm = TRUE
      )
    },
    integer(1)
  )

  # ------------------------------------------------------------------
  # Prediction template:
  # one row for every time value x treatment group combination
  # ------------------------------------------------------------------

  prediction_template <- data.table::data.table(
    id = seq_len(2L * length(month_values))
  )

  prediction_template[[timeVar]] <- rep(
    month_values,
    times = 2L
  )

  prediction_template[[group_col]] <- factor(
    rep(
      c(treated_group, untreated_group),
      each = length(month_values)
    ),
    levels = levels(df[[group_col]])
  )

  # Preserve factor levels for the time variable if it is a factor
  if (is.factor(df[[timeVar]])) {
    prediction_template[[timeVar]] <- factor(
      prediction_template[[timeVar]],
      levels = levels(df[[timeVar]])
    )
  }

  # ------------------------------------------------------------------
  # Prediction on LINK scale
  #
  # For a log-binomial model:
  #
  #     eta = log(p)
  #
  # ------------------------------------------------------------------

  eta <- stats::predict(
    fit,
    newdata = prediction_template,
    type = "link"
  )

  eta <- as.numeric(eta)

  # ------------------------------------------------------------------
  # Design matrix for new observations
  # ------------------------------------------------------------------

  Xnew <- stats::model.matrix(
    stats::delete.response(stats::terms(fit)),
    prediction_template,
    contrasts.arg = fit$contrasts
  )

  # Make sure columns have exactly the same order as coefficients
  coef_names <- names(stats::coef(fit))

  missing_coef <- setdiff(coef_names, colnames(Xnew))

  if (length(missing_coef) > 0L) {
    stop(
      paste0(
        "Prediction design matrix is missing coefficient columns: ",
        paste(missing_coef, collapse = ", ")
      )
    )
  }

  Xnew <- Xnew[, coef_names, drop = FALSE]

  # ------------------------------------------------------------------
  # Robust covariance matrix of GEE coefficients
  # ------------------------------------------------------------------

  V <- stats::vcov(fit)

  V <- V[
    coef_names,
    coef_names,
    drop = FALSE
  ]

  # ------------------------------------------------------------------
  # Standard error on LINK scale:
  #
  # SE(eta) = sqrt(x' V x)
  # ------------------------------------------------------------------

  se_eta <- sqrt(
    rowSums(
      (Xnew %*% V) * Xnew
    )
  )

  # ------------------------------------------------------------------
  # Convert from log scale to probability scale
  #
  # log(p) = eta
  # p      = exp(eta)
  # ------------------------------------------------------------------

  predicted_prob <- exp(eta)

  # ------------------------------------------------------------------
  # Delta-method SE on probability scale
  #
  # If:
  #
  #     p = exp(eta)
  #
  # then:
  #
  #     dp/deta = exp(eta) = p
  #
  # therefore:
  #
  #     SE(p) ≈ p * SE(eta)
  # ------------------------------------------------------------------

  prediction_se <- predicted_prob * se_eta

  # ------------------------------------------------------------------
  # Confidence interval
  #
  # Construct CI on log/link scale first, then transform.
  # ------------------------------------------------------------------

  z_value <- stats::qnorm(
    1 - (1 - conf_level) / 2
  )

  lower_eta <- eta - z_value * se_eta
  upper_eta <- eta + z_value * se_eta

  ci_lb <- exp(lower_eta)
  ci_ub <- exp(upper_eta)

  # ------------------------------------------------------------------
  # Store prediction components
  # ------------------------------------------------------------------
  prediction_summary <- data.table::copy(prediction_template)
  prediction_summary[, `:=`(
    n = c(treated_counts, untreated_counts),
    mean = predicted_prob,
    se = prediction_se,
    ci_lb = ci_lb,
    ci_ub = ci_ub
  )]

  return(prediction_summary[])
}