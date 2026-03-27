# --- helper wrapper function to fit a cox model, return messages to log
fitmod_cox <- function(model_formula, iptw = "ip_weights",
                       model_type = "adj", aesi_name = NULL, aesifup_input, ...) {
  tryCatch(
    {
      if (model_type == "crude") {
        return(survival::coxph(model_formula,
          cluster = person_id, id = person_id,
          robust = TRUE, data = aesifup_input
        ))
      }
      if (model_type == "adj") {
        w <- aesifup_input[[iptw]]
        pid <- aesifup_input[["person_id"]]

        return(survival::coxph(
          model_formula,
          data = aesifup_input,
          robust = TRUE,
          id = pid,
          cluster = pid,
          weights = w
        ))
      }
    },
    error = function(cond) {
      logger::log_info(paste("Cox model fitting returning error for ", model_type, " model, ", aesi_name))
      logger::log_info(conditionMessage(cond))
    },
    warning = function(cond) {
      logger::log_info(paste("Cox model fitting returning warning for ", model_type, " model, ", aesi_name))
      logger::log_info(conditionMessage(cond))
    }
  )
}
