#' Fit a log-binomial regression model for risk ratio estimation
#'
#' Fits a generalised linear model with a binomial family and log link to
#' directly estimate risk ratios (and their Wald confidence intervals) in
#' matched or weighted cohort data.  Uses only base-R `glm()` — no
#' additional package required.
#'
#' @param model_formula formula passed to `glm`, e.g. `eventCount ~ group`
#' @param model_type `"crude"` (unweighted) or `"adj"` (IPTW-weighted)
#' @param iptw column name of the inverse-probability weight variable
#' @param aesi_name outcome label used in log messages
#' @param aesifup_input data.frame / data.table containing the analysis dataset
#' @param dummy_code scalar returned when model fitting fails
#'
#' @return a `glm` object on success, or a character error/warning message
#'   on failure (handled gracefully downstream by `extract_ests_logbin`)
#' @export
fitmod_logbin <- function(model_formula,
                          model_type  = NULL,
                          iptw        = "ip_weight",
                          aesi_name   = NULL,
                          aesifup_input,
                          dummy_code  = -88) {
  # Use withCallingHandlers for warnings so glm() can complete even when it
  # emits a convergence warning (common with binomial/log link + weights).
  # tryCatch is reserved for hard errors only.
  tryCatch({
    withCallingHandlers({
      df <- as.data.frame(aesifup_input)
      if (!(model_type %in% c("crude", "adj"))) {
        stop("model_type must be either 'crude' or 'adj'")
      }

      glm_args <- list(
        formula = model_formula,
        family = binomial(link = "log"),
        data = df
      )

      if (model_type == "adj") {
        if (is.null(iptw) || !(iptw %in% names(df))) {
          stop(paste0("Weight column '", iptw, "' not found in aesifup_input"))
        }
        # Attach weights as a column so glm() can find it by name in data
        df[[".wt"]] <- df[[iptw]]
        glm_args$data <- df
        glm_args$weights <- df[[".wt"]]
      }

      do.call(stats::glm, glm_args)
    },
    warning = function(w) {
      logger::log_info(paste("Log-binomial warning for", model_type,
                             "model,", aesi_name, ":", conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
  },
  error = function(cond) {
    msg <- paste("Log-binomial model returning error for", model_type,
                 "model,", aesi_name)
    logger::log_info(msg)
    logger::log_info(conditionMessage(cond))
    return(msg)
  })
}
