#' Fit Poisson fallback for time-varying GEE
#'
#' Re-fits a time-varying comparative model with a Poisson log link when the
#' primary binomial-log GEE fit fails, while preserving the same clustered GEE
#' structure used by the main model.
#'
#' @param model_formula_input Model formula to fit.
#' @param model_type Either `"crude"` or `"adj"`.
#' @param model_data Input data for the fit.
#' @param idCol Name of the clustering id column.
#' @param iptw Weight column name for adjusted fits.
#' @param target_aesi Outcome label used for logging.
#' @param binomial_error Original binomial model error message.
#'
#' @return A fitted `geeglm` object or a character error message.
#' @keywords internal
.fit_timevarying_poisson <- function(model_formula_input,
                                     model_type,
                                     model_data,
                                     idCol = "person_id_num",
                                     iptw = NULL,
                                     target_aesi = NULL,
                                     binomial_error = NULL) {
  tryCatch({
    withCallingHandlers({
      df <- as.data.frame(model_data)

      if (!(idCol %in% names(df))) {
        stop(paste0("Column '", idCol, "' not found in aesifup_input"))
      }

      weight_vec <- NULL
      if (model_type == "adj") {
        if (is.null(iptw) || !(iptw %in% names(df))) {
          stop(paste0("Weight column '", iptw, "' not found in aesifup_input"))
        }
        weight_vec <- df[[iptw]]
      }

      geeglm_args <- list(
        formula = model_formula_input,
        data = df,
        family = stats::poisson(link = "log"),
        id = df[[idCol]]
      )
      if (!is.null(weight_vec)) {
        geeglm_args$weights <- weight_vec
      }

      poisson_fit <- do.call(geepack::geeglm, geeglm_args)
      attr(poisson_fit, "timevarying_family") <- "poisson"
      attr(poisson_fit, "timevarying_fallback") <- TRUE
      attr(poisson_fit, "binomial_error") <- binomial_error
      poisson_fit
    },
    warning = function(cond) {
      warning_message <- paste("Time-varying Poisson-log GEE warning for",
                               model_type, "model,", target_aesi)
      logger::log_info(warning_message)
      logger::log_info(conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
  },
  error = function(cond) {
    error_message <- paste("Time-varying Poisson-log GEE error for",
                           model_type, "model,", target_aesi)
    logger::log_info(error_message)
    logger::log_info(conditionMessage(cond))
    return(error_message)
  })
}