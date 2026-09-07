# internal function to fit time-varying binomial GEE models
# uses geepack::geeglm with log link and clustered id
fitmod_gee_tv <- function(model_formula,
                          idCol = "person_id_num",
                          iptw = "ip_weight",
                          model_type = NULL,
                          aesi_name = NULL,
                          aesifup_input,
                          ...){
  tryCatch({
    withCallingHandlers({
      df <- as.data.frame(aesifup_input)

      if (!(idCol %in% names(df))) {
        stop(paste0("Column '", idCol, "' not found in aesifup_input"))
      }

      if (!(model_type %in% c("crude", "adj"))) {
        stop("model_type must be either 'crude' or 'adj'")
      }

      weight_vec <- NULL
      if (model_type == "adj") {
        if (is.null(iptw) || !(iptw %in% names(df))) {
          stop(paste0("Weight column '", iptw, "' not found in aesifup_input"))
        }
        weight_vec <- df[[iptw]]
      }

      geeglm_args <- list(
        formula = model_formula,
        data = df,
        family = binomial(link = "log"),
        id = df[[idCol]]
      )
      if (!is.null(weight_vec)) {
        geeglm_args$weights <- weight_vec
      }
      return(do.call(geepack::geeglm, geeglm_args))
    },
    warning = function(cond){
      warning_message <- paste("Time-varying GEE warning for", model_type,
                               "model,", aesi_name)
      logger::log_info(warning_message)
      logger::log_info(conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
  },
  error = function(cond){
    error_message <- paste("Time-varying GEE error for", model_type,
                           "model,", aesi_name)
    logger::log_info(error_message)
    logger::log_info(conditionMessage(cond))
    return(error_message)
  })
}