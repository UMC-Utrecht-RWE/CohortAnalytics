#' Build an empty time-varying prediction table
#'
#' Creates an empty result table for time-varying prediction summaries.
#'
#' @param timeVar Name of the time index column.
#'
#' @return An empty data.table with the expected treated and untreated
#'   prediction columns.
#' @keywords internal
.empty_timevarying_prediction_table <- function(timeVar = "time_value") {
  out <- data.table::data.table(
    time_value = numeric(),
    treated_n = integer(),
    treated_mean = numeric(),
    treated_sd = numeric(),
    treated_ci_lb = numeric(),
    treated_ci_ub = numeric(),
    untreated_n = integer(),
    untreated_mean = numeric(),
    untreated_sd = numeric(),
    untreated_ci_lb = numeric(),
    untreated_ci_ub = numeric()
  )

  data.table::setnames(out, "time_value", timeVar)
  out
}