# extract comparative estimates from a time-varying binomial GEE model
# returns RR and an approximate RD at the reference time level
extract_ests_gee_tv <- function(modelobj,
                                scale_IR = 100000,
                                throw_warning = FALSE) {
  if (is.character(modelobj)) {
    if (throw_warning) {
      ("warning: model object is character string, returning NA flags -88")
    }
    return(data.frame(irr_est = -88,
                      irr_lb = -88,
                      irr_ub = -88,
                      ird_est = -88,
                      ird_lb = -88,
                      ird_ub = -88))
  }

  coef_vec <- tryCatch(stats::coef(modelobj), error = function(e) NULL)
  if (is.null(coef_vec)) {
    return(data.frame(irr_est = -88,
                      irr_lb = -88,
                      irr_ub = -88,
                      ird_est = -88,
                      ird_lb = -88,
                      ird_ub = -88))
  }

  coef_names <- names(coef_vec)
  int_idx <- which(coef_names == "(Intercept)")
  group_idx <- grep("^group.*EXPOSED$|EXPOSED", coef_names)

  if (length(int_idx) == 0 || length(group_idx) == 0) {
    return(data.frame(irr_est = -88,
                      irr_lb = -88,
                      irr_ub = -88,
                      ird_est = -88,
                      ird_lb = -88,
                      ird_ub = -88))
  }

  int_idx <- int_idx[1]
  group_idx <- group_idx[1]

  vcov_mat <- tryCatch(stats::vcov(modelobj), error = function(e) NULL)
  if (is.null(vcov_mat) && !is.null(modelobj$var)) {
    vcov_mat <- modelobj$var
  }
  if (is.null(vcov_mat)) {
    return(data.frame(irr_est = -88,
                      irr_lb = -88,
                      irr_ub = -88,
                      ird_est = -88,
                      ird_lb = -88,
                      ird_ub = -88))
  }

  log_rr <- coef_vec[group_idx]
  log_rr_se <- sqrt(vcov_mat[group_idx, group_idx])

  rr_out <- data.frame(
    irr_est = exp(log_rr),
    irr_lb = exp(log_rr + qnorm(0.025) * log_rr_se),
    irr_ub = exp(log_rr + qnorm(0.975) * log_rr_se)
  )

  mean2 <- c(coef_vec[int_idx], coef_vec[group_idx])
  vcov2 <- vcov_mat[c(int_idx, group_idx), c(int_idx, group_idx), drop = FALSE]

  rd_est <- exp(mean2[1] + mean2[2]) - exp(mean2[1])
  rd_se <- msm::deltamethod(~exp(x1 + x2) - exp(x1), mean2, vcov2)

  rd_out <- data.frame(
    ird_est = rd_est,
    ird_lb = rd_est + qnorm(0.025) * rd_se,
    ird_ub = rd_est + qnorm(0.975) * rd_se
  ) * scale_IR

  outobj <- data.frame(rr_out, rd_out)
  rownames(outobj) <- NULL
  outobj
}


  # Replace the ratio portion returned by extract_ests_gee_tv with the
  # coefficient-based RR and robust Wald CI. This works for both the primary
  # binomial-log GEE and the Poisson-log fallback.
  .set_timevarying_rr <- function(extracted_estimates,
                                  model_object,
                                  conf_level = 0.95) {
    if (!inherits(model_object, "geeglm")) {
      return(extracted_estimates)
    }

    coef_table <- summary(model_object)$coefficients
    group_rows <- grep("^group", rownames(coef_table))
    if (length(group_rows) != 1L) {
      logger::log_info(paste(
        "Could not identify one group coefficient in time-varying GEE for",
        target_aesi
      ))
      return(extracted_estimates)
    }

    estimate_col <- if ("Estimate" %in% colnames(coef_table)) {
      "Estimate"
    } else {
      1L
    }
    se_col <- if ("Std.err" %in% colnames(coef_table)) {
      "Std.err"
    } else if ("Std. Error" %in% colnames(coef_table)) {
      "Std. Error"
    } else {
      2L
    }

    beta_group <- coef_table[group_rows, estimate_col]
    robust_se <- coef_table[group_rows, se_col]
    z_value <- stats::qnorm(1 - (1 - conf_level) / 2)

    rr_values <- c(
      exp(beta_group),
      exp(beta_group - z_value * robust_se),
      exp(beta_group + z_value * robust_se)
    )

    # The first three columns from extract_ests_gee_tv are the ratio estimate,
    # lower confidence limit, and upper confidence limit.
    extracted_estimates[1, 1:3] <- rr_values
    extracted_estimates
  }