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
