# internal function to wrangle model fitting for poisson outcomes
# supports geeM::geem() (default) with GEE and biglm::bigglm() as memory-efficient fallback
# if error or warning, returns string
# relies on availability of iptw column, variable names, person_id_num, aesifup in environment, etc.
fitmod_gee <- function(model_formula, iptw = "iptw", model_type = NULL, aesi_name = NULL, aesifup_input, model_library = "geeM", chunksize = 100000, ...) {
  logger::log_info(paste0("Fitting model with ", model_library, " for ", model_type, " model, ", aesi_name))

  # Issue notice upfront if using bigglm fallback
  if (model_library == "bigglm") {
    logger::log_warn(paste0(
      "Using bigglm with Poisson family for ", aesi_name,
      " - Note: Poisson GLM with sandwich robust SEs is used instead of GEE Poisson; clustering structure is not explicitly modeled"
    ))
  }

  tryCatch(
    {
      if (model_library == "geeM") {
        logger::log_info(paste0("Using geeM::geem() for ", aesi_name, " with model type: ", model_type))
        if (model_type == "crude") {
          return(geeM::geem(
            formula = model_formula, data = aesifup_input, family = "poisson",
            id = person_id_num
          ))
        }
        if (model_type == "adj") {
          return(geeM::geem(
            formula = model_formula, data = aesifup_input, family = "poisson",
            id = person_id_num, weights = aesifup_input[[iptw]]
          ))
        }
      } else if (model_library == "bigglm") {
        logger::log_info(paste0("Using biglm::bigglm() for ", aesi_name, " with model type: ", model_type))
        # Use biglm with sandwich robust SEs for memory efficiency on large datasets
        # Weights must be specified as one-sided formula per biglm documentation
        weight_formula <- if (model_type == "adj") as.formula(paste0("~", iptw)) else NULL

        if (model_type == "crude") {
          return(biglm::bigglm(
            formula = model_formula,
            data = aesifup_input,
            family = poisson(),
            chunksize = chunksize,
            maxit = 100,
            tolerance = 1e-6,
            sandwich = TRUE
          ))
        }
        if (model_type == "adj") {
          return(biglm::bigglm(
            formula = model_formula,
            data = aesifup_input,
            weights = weight_formula,
            family = poisson(),
            chunksize = chunksize,
            maxit = 100,
            tolerance = 1e-6,
            sandwich = TRUE
          ))
        }
      } else {
        stop(paste0("Unknown model_library: ", model_library, ". Supported options are 'geeM' and 'bigglm'."))
      }
    },
    error = function(cond) {
      error_message <- paste("Model fitting returning error for ", model_type, " model, ", aesi_name)
      logger::log_info(error_message)
      logger::log_info(conditionMessage(cond))
      return(error_message)
    },
    warning = function(cond) {
      warning_message <- paste("Model fitting returning warning for ", model_type, " model, ", aesi_name)
      logger::log_info(warning_message)
      logger::log_info(conditionMessage(cond))
      return(warning_message)
    }
  )
}
