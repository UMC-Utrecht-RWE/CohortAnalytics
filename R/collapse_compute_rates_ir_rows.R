#' Collapse time-varying rows for incidence summaries
#'
#' Reduces repeated person-period rows to one row per person for descriptive
#' incidence counts in the time-varying incidence path. If a person has an event,
#' the first event row is retained; otherwise the final observed row is kept.
#'
#' @param group_rows Data frame or data.table containing one treatment group.
#' @param type Analysis type passed from `compute_rates_cohort`.
#' @param incidence_model Incidence model selector.
#' @param timeVar Name of the time index column.
#' @param idCol Name of the person identifier column.
#' @param eventCol Name of the event count column.
#'
#' @return A data.frame-like object with at most one row per person.
#' @keywords internal
.collapse_compute_rates_ir_rows <- function(group_rows,
                                            type,
                                            incidence_model,
                                            timeVar,
                                            idCol,
                                            eventCol) {
  use_timevarying_ir_rows <- (
    type == "incidence" &&
      incidence_model == "timevarying" &&
      !is.null(timeVar) &&
      idCol %in% names(group_rows) &&
      timeVar %in% names(group_rows)
  )

  if (!use_timevarying_ir_rows) {
    return(group_rows)
  }

  row_order <- order(group_rows[[idCol]], group_rows[[timeVar]], na.last = TRUE)
  group_rows <- group_rows[row_order, , drop = FALSE]
  row_ids <- split(seq_len(nrow(group_rows)), group_rows[[idCol]], drop = TRUE)

  keep_rows <- vapply(row_ids,
                      FUN.VALUE = integer(1),
                      function(person_rows) {
                        event_values <- group_rows[[eventCol]][person_rows]
                        event_hits <- which(!is.na(event_values) & event_values > 0)

                        if (length(event_hits) > 0) {
                          person_rows[event_hits[1]]
                        } else {
                          utils::tail(person_rows, 1)
                        }
                      })

  group_rows[keep_rows, , drop = FALSE]
}