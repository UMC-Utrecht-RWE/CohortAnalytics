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

      if (model_type == "crude") {
        return(geepack::geeglm(
          formula = model_formula,
          data = df,
          family = binomial(link = "log"),
          id = df[[idCol]]
        ))
      }

      if (model_type == "adj") {
        if (!(iptw %in% names(df))) {
          stop(paste0("Weight column '", iptw, "' not found in aesifup_input"))
        }
        return(geepack::geeglm(
          formula = model_formula,
          data = df,
          family = binomial(link = "log"),
          id = df[[idCol]],
          weights = df[[iptw]]
        ))
      }

      stop("model_type must be either 'crude' or 'adj'")
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
