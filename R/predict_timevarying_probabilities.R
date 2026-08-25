#' Predict counterfactual time-varying probabilities
#'
#' Creates treated and untreated copies of the model dataset, subsets each copy
#' to every observed time value, and returns summary statistics for the
#' corresponding predicted probabilities from a fitted time-varying model.
#'
#' @param fit A fitted time-varying model object.
#' @param model_data Data used as the prediction template.
#' @param timeVar Name of the time index column.
#' @param idCol Name of the person identifier column.
#' @param group_col Name of the treatment group column.
#' @param treated_group Value representing the treated strategy.
#' @param untreated_group Value representing the untreated strategy.
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
    return(data.table::data.table(
      time_value = numeric(),
      treated_n = integer(),
      treated_mean = numeric(),
      treated_ci_lb = numeric(),
      treated_ci_ub = numeric(),
      untreated_n = integer(),
      untreated_mean = numeric(),
      untreated_ci_lb = numeric(),
      untreated_ci_ub = numeric()
    ))
  }

  df <- as.data.frame(model_data)

  if (is.null(timeVar) || !(timeVar %in% names(df))) {
    stop("A valid timeVar column is required for time-varying predictions")
  }

  if (!(group_col %in% names(df))) {
    stop(paste0("Column '", group_col, "' not found in model_data"))
  }

  dat1 <- df
  dat0 <- df

  if (is.factor(df[[group_col]])) {
    group_levels <- levels(df[[group_col]])
    dat1[[group_col]] <- factor(rep(treated_group, nrow(dat1)), levels = group_levels)
    dat0[[group_col]] <- factor(rep(untreated_group, nrow(dat0)), levels = group_levels)
  } else {
    dat1[[group_col]] <- treated_group
    dat0[[group_col]] <- untreated_group
  }

  month_values <- unique(df[[timeVar]][!is.na(df[[timeVar]])])
  month_predictions <- vector("list", length(month_values))

  for (month_index in seq_along(month_values)) {
    month_value <- month_values[[month_index]]
    keep_rows <- !is.na(df[[timeVar]]) & df[[timeVar]] == month_value

    month_dat1 <- dat1[keep_rows, , drop = FALSE]
    month_dat0 <- dat0[keep_rows, , drop = FALSE]

    p1 <- tryCatch(stats::predict(fit, newdata = month_dat1, type = "response"),
                   error = function(cond) cond)
    p0 <- tryCatch(stats::predict(fit, newdata = month_dat0, type = "response"),
                   error = function(cond) cond)

    p1_summary <- if (inherits(p1, "error")) NULL else .summarize_timevarying_predictions(p1, conf_level = conf_level)
    p0_summary <- if (inherits(p0, "error")) NULL else .summarize_timevarying_predictions(p0, conf_level = conf_level)

    month_predictions[[month_index]] <- data.table::data.table(
      time_value = month_value,
      treated_n = if (is.null(p1_summary)) 0L else p1_summary$n,
      treated_mean = if (is.null(p1_summary)) NA_real_ else p1_summary$mean,
      treated_ci_lb = if (is.null(p1_summary)) NA_real_ else p1_summary$ci_lb,
      treated_ci_ub = if (is.null(p1_summary)) NA_real_ else p1_summary$ci_ub,
      untreated_n = if (is.null(p0_summary)) 0L else p0_summary$n,
      untreated_mean = if (is.null(p0_summary)) NA_real_ else p0_summary$mean,
      untreated_ci_lb = if (is.null(p0_summary)) NA_real_ else p0_summary$ci_lb,
      untreated_ci_ub = if (is.null(p0_summary)) NA_real_ else p0_summary$ci_ub
    )
  }

  prediction_summary <- data.table::rbindlist(month_predictions, use.names = TRUE, fill = TRUE)
  data.table::setnames(prediction_summary, "time_value", timeVar)
  prediction_summary
}