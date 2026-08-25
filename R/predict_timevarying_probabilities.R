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
.predict_timevarying_probabilities <- function(fit,
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

  treated_counts <- vapply(
    month_values,
    function(month_value) {
      sum(df[[timeVar]] == month_value & df[[group_col]] == treated_group, na.rm = TRUE)
    },
    integer(1)
  )
  untreated_counts <- vapply(
    month_values,
    function(month_value) {
      sum(df[[timeVar]] == month_value & df[[group_col]] == untreated_group, na.rm = TRUE)
    },
    integer(1)
  )

  prediction_template <- df[rep(1, length(month_values) * 2L), , drop = FALSE]
  prediction_template[[timeVar]] <- rep(month_values, times = 2L)

  if (is.factor(df[[group_col]])) {
    prediction_template[[group_col]] <- factor(
      rep(c(treated_group, untreated_group), each = length(month_values)),
      levels = levels(df[[group_col]])
    )
  } else {
    prediction_template[[group_col]] <- rep(c(treated_group, untreated_group), each = length(month_values))
  }

  predicted_values <- tryCatch(
    stats::predict(fit, newdata = prediction_template, type = "response"),
    error = function(cond) cond
  )

  if (inherits(predicted_values, "error")) {
    return(.empty_timevarying_prediction_table(timeVar = timeVar))
  }

  link_values <- tryCatch(
    stats::predict(fit, newdata = prediction_template, type = "link"),
    error = function(cond) cond
  )

  vcov_mat <- tryCatch(stats::vcov(fit), error = function(cond) NULL)
  if (is.null(vcov_mat) && !is.null(fit$var)) {
    vcov_mat <- fit$var
  }

  se_link <- rep(NA_real_, length(predicted_values))
  z_value <- stats::qnorm(1 - (1 - conf_level) / 2)

  if (!inherits(link_values, "error") && !is.null(vcov_mat)) {
    model_terms <- stats::delete.response(stats::terms(fit))
    model_matrix <- stats::model.matrix(object = model_terms,
                                        data = prediction_template,
                                        contrasts.arg = fit$contrasts)

    vcov_mat <- vcov_mat[colnames(model_matrix), colnames(model_matrix), drop = FALSE]
    se_link <- sqrt(pmax(rowSums((model_matrix %*% vcov_mat) * model_matrix), 0))
  }

  prediction_se <- predicted_values * se_link
  ci_lb <- if (inherits(link_values, "error")) {
    rep(NA_real_, length(predicted_values))
  } else {
    exp(link_values - z_value * se_link)
  }
  ci_ub <- if (inherits(link_values, "error")) {
    rep(NA_real_, length(predicted_values))
  } else {
    exp(link_values + z_value * se_link)
  }

  month_count <- length(month_values)
  treated_predictions <- predicted_values[seq_len(month_count)]
  untreated_predictions <- predicted_values[month_count + seq_len(month_count)]
  treated_se <- prediction_se[seq_len(month_count)]
  untreated_se <- prediction_se[month_count + seq_len(month_count)]
  treated_ci_lb <- ci_lb[seq_len(month_count)]
  treated_ci_ub <- ci_ub[seq_len(month_count)]
  untreated_ci_lb <- ci_lb[month_count + seq_len(month_count)]
  untreated_ci_ub <- ci_ub[month_count + seq_len(month_count)]

  prediction_summary <- data.table::data.table(
    time_value = month_values,
    treated_n = treated_counts,
    treated_mean = treated_predictions,
    treated_sd = treated_se,
    treated_ci_lb = treated_ci_lb,
    treated_ci_ub = treated_ci_ub,
    untreated_n = untreated_counts,
    untreated_mean = untreated_predictions,
    untreated_sd = untreated_se,
    untreated_ci_lb = untreated_ci_lb,
    untreated_ci_ub = untreated_ci_ub
  )

  data.table::setnames(prediction_summary, "time_value", timeVar)
  prediction_summary
}