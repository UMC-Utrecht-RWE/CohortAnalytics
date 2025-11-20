# internal function to wrangle geeM:geem() model fiting
# fits the conditional poisson model; if error or warning, returns string
# relies on availability of iptw column, variable names, person_id_num, aesifup in environment, etc.
fitmod_gee <- function(model_formula, iptw = "iptw", model_type = NULL,aesi_name = NULL,aesifup_input, ...){
  tryCatch({
    if(model_type == "crude"){
      return(geeM::geem(formula = model_formula, data = aesifup_input, family = "poisson",
                        id = person_id_num))
    }
    if(model_type == "adj"){
      return(geeM::geem(formula = model_formula, data = aesifup_input, family = "poisson",
                        id = person_id_num, weights = aesifup_input[[iptw]]))
    }
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
