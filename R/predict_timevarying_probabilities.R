#' Predict counterfactual time-varying probabilities
#'
#' Creates treated and untreated copies of the model dataset, subsets each copy
#' to every observed time value, and returns the corresponding predicted
#' probabilities from a fitted time-varying model.
#'
#' @param fit A fitted time-varying model object.
#' @param model_data Data used as the prediction template.
#' @param timeVar Name of the time index column.
#' @param idCol Name of the person identifier column.
#' @param group_col Name of the treatment group column.
#' @param treated_group Value representing the treated strategy.
#' @param untreated_group Value representing the untreated strategy.
#'
#' @return A list containing the time variable name and per-timepoint predicted
#'   probabilities under treated and untreated assignments.
#' @keywords internal
.predict_timevarying_probabilities <- function(fit,
                                               model_data,
                                               timeVar,
                                               idCol = "person_id_num",
                                               group_col = "group",
                                               treated_group = "EXPOSED",
                                               untreated_group = "CONTROL") {
  if (is.character(fit)) {
    return(list(error = fit,
                time_variable = timeVar,
                months = list()))
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
  month_names <- as.character(month_values)
  month_predictions <- setNames(vector("list", length(month_values)), month_names)

  for (month_index in seq_along(month_values)) {
    month_value <- month_values[[month_index]]
    keep_rows <- !is.na(df[[timeVar]]) & df[[timeVar]] == month_value

    month_dat1 <- dat1[keep_rows, , drop = FALSE]
    month_dat0 <- dat0[keep_rows, , drop = FALSE]

    p1 <- tryCatch(stats::predict(fit, newdata = month_dat1, type = "response"),
                   error = function(cond) cond)
    p0 <- tryCatch(stats::predict(fit, newdata = month_dat0, type = "response"),
                   error = function(cond) cond)

    month_predictions[[month_index]] <- list(
      month_value = month_value,
      row_index = which(keep_rows),
      person_id = if (idCol %in% names(df)) df[[idCol]][keep_rows] else NULL,
      p_exposed = if (inherits(p1, "error")) NULL else p1,
      p_control = if (inherits(p0, "error")) NULL else p0,
      mean_p_exposed = if (inherits(p1, "error")) NA_real_ else mean(p1, na.rm = TRUE),
      mean_p_control = if (inherits(p0, "error")) NA_real_ else mean(p0, na.rm = TRUE),
      error = c(if (inherits(p1, "error")) conditionMessage(p1) else NULL,
                if (inherits(p0, "error")) conditionMessage(p0) else NULL)
    )
  }

  list(time_variable = timeVar,
       months = month_predictions)
}