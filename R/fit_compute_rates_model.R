#' Fit one compute-rates model
#'
#' Routes `compute_rates_cohort` model fitting to the standard GEE, log-binomial,
#' or time-varying GEE implementation and applies the Poisson fallback for the
#' time-varying path when the binomial-log fit fails.
#'
#' @param model_formula_input Model formula to fit.
#' @param model_type Either `"crude"` or `"adj"`.
#' @param model_data Input data for the fit.
#' @param use_timevarying Logical flag for the time-varying GEE path.
#' @param use_logbin Logical flag for the log-binomial path.
#' @param iptw Weight column name for adjusted fits.
#' @param idCol Name of the clustering id column.
#' @param target_aesi Outcome label used for logging.
#'
#' @return A fitted model object or a character error message, matching the
#'   downstream behavior expected by `compute_rates_cohort`.
#' @keywords internal
.fit_compute_rates_model <- function(model_formula_input,
                                     model_type,
                                     model_data,
                                     use_timevarying,
                                     use_logbin,
                                     iptw,
                                     idCol,
                                     target_aesi) {
  if (!(model_type %in% c("crude", "adj"))) {
    stop("model_type must be either 'crude' or 'adj'")
  }

  iptw_arg <- if (model_type == "adj") iptw else NULL

  if (use_timevarying) {
    tv_model_data <- as.data.frame(model_data)

    primary_fit <- fitmod_gee_tv(model_formula_input,
                                 model_type = model_type,
                                 iptw = iptw_arg,
                                 idCol = idCol,
                                 aesi_name = target_aesi,
                                 aesifup_input = tv_model_data)

    if (is.character(primary_fit) && length(primary_fit) == 1L) {
      deviation_message <- paste(
        "Time-varying GEE model deviation for", model_type, "model,",
        target_aesi, ": binomial-log failed; using Poisson-log fallback"
      )
      logger::log_info(deviation_message)
      logger::log_info(primary_fit)

      return(.fit_timevarying_poisson(
        model_formula_input = model_formula_input,
        model_type = model_type,
        model_data = tv_model_data,
        idCol = idCol,
        iptw = iptw_arg,
        target_aesi = target_aesi,
        binomial_error = primary_fit
      ))
    }

    if (inherits(primary_fit, "geeglm")) {
      attr(primary_fit, "timevarying_family") <- "binomial"
      attr(primary_fit, "timevarying_fallback") <- FALSE
    }
    return(primary_fit)
  }

  if (use_logbin) {
    return(fitmod_logbin(model_formula_input,
                         model_type = model_type,
                         iptw = iptw_arg,
                         aesi_name = target_aesi,
                         aesifup_input = model_data))
  }

  fitmod_gee(model_formula_input,
             model_type = model_type,
             iptw = iptw_arg,
             idCol = idCol,
             aesi_name = target_aesi,
             aesifup_input = model_data)
}