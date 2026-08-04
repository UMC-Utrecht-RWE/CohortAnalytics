# internal function to wrangle geeM:geem() model fiting
# fits the conditional poisson model; if error or warning, returns string
# relies on availability of iptw column, variable names, person_id_num, aesifup in environment, etc.
fitmod_gee <- function(model_formula, iptw = "iptw", model_type = NULL,
                       idCol = "person_id_num", aesi_name = NULL, aesifup_input, ...){
  tryCatch({
    if (!(model_type %in% c("crude", "adj"))) {
      stop("model_type must be either 'crude' or 'adj'")
    }

    if (!(idCol %in% names(aesifup_input))) {
      stop(paste0("Column '", idCol, "' not found in aesifup_input"))
    }

    id_vec <- aesifup_input[[idCol]]
    weight_vec <- NULL
    if (model_type == "adj") {
      if (is.null(iptw) || !(iptw %in% names(aesifup_input))) {
        stop(paste0("Weight column '", iptw, "' not found in aesifup_input"))
      }
      weight_vec <- aesifup_input[[iptw]]
    }

    geem_args <- list(
      formula = model_formula,
      data = aesifup_input,
      family = "poisson",
      id = id_vec
    )
    if (!is.null(weight_vec)) {
      geem_args$weights <- weight_vec
    }
    return(do.call(geeM::geem, geem_args))
  },
  error = function(cond){
    error_message <- paste("Model fitting returning error for ",model_type, " model, ", aesi_name)
    logger::log_info(error_message)
    logger::log_info(conditionMessage(cond))
    return(error_message)
  },
  warning = function(cond){
    warning_message <- paste("Model fitting returning warning for ",model_type, " model, ", aesi_name)
    logger::log_info(warning_message)
    logger::log_info(conditionMessage(cond))
    return(warning_message)
  }
  )
}
