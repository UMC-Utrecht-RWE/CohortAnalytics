#' Fit one-group compute-rates model
#'
#' Fits the exposed-only or control-only model used for model-based rate
#' estimation in `compute_rates_cohort`.
#'
#' @param group_name Group to subset to before fitting.
#' @param model_type Either `"crude"` or `"adj"`.
#' @param aesifup_input Full analysis input dataset.
#' @param use_timevarying Logical flag for the time-varying GEE path.
#' @param model_formula_no_group Formula used for non-time-varying one-group fits.
#' @param model_formula_no_group_tv Optional formula used for time-varying one-group fits.
#' @param use_logbin Logical flag for the log-binomial path.
#' @param iptw Weight column name for adjusted fits.
#' @param idCol Name of the clustering id column.
#' @param target_aesi Outcome label used for logging.
#'
#' @return A fitted model object or a character error message.
#' @keywords internal
.fit_compute_rates_group_model <- function(group_name,
                                           model_type,
                                           aesifup_input,
                                           use_timevarying,
                                           model_formula_no_group,
                                           model_formula_no_group_tv = NULL,
                                           use_logbin,
                                           iptw,
                                           idCol,
                                           target_aesi) {
  group_data <- aesifup_input[aesifup_input$group == group_name, ]
  model_formula_group <- if (use_timevarying) {
    model_formula_no_group_tv
  } else {
    model_formula_no_group
  }

  .fit_compute_rates_model(model_formula_input = model_formula_group,
                           model_type = model_type,
                           model_data = group_data,
                           use_timevarying = use_timevarying,
                           use_logbin = use_logbin,
                           iptw = iptw,
                           idCol = idCol,
                           target_aesi = target_aesi)
}