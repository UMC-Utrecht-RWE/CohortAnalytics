#' Wrapper function to compute statistics relating to prevelance
#'
#'  rates, ratios and differences
#'  per group and between groups
#'  with or without weights
#'
#' key functionality:
#'  est_inc_prev
#'  create_risk_table
#'  estimate_HR
#'
#' output:
#'  point estimates (mean)
#'  CI based on bootstrap results
#'
#' @param aesifup_input data
#' @param eventCol column of event 0/1
#' @param groupCol column identifying exposed/control group
#' @param iptw column of weight
#' @param bootstrap bootstrap result to calculate CI
#' @param adjust either 'adj' adjusted or 'unadj' unadjusted by iptw
#' @param estimate_survival TRUE/FALSE whether to apply survival model
#' @param risk_window risk window for survival analysis
#'
#' @export
compute_prev_rates <- function(
    aesifup_input,
    eventCol = "event",
    groupCol = "group",
    iptw = "wt",
    bootstrap,
    adjust = "adj",
    scale_IR = 1,
    estimate_survival = FALSE,
    risk_window = NULL
) {
  n_pat_exp <- aesifup_input[aesifup_input[[groupCol]] == "EXPOSED", .N]
  n_pat_con <- aesifup_input[aesifup_input[[groupCol]] == "CONTROL", .N]
  n_out_exp <- aesifup_input[aesifup_input[[groupCol]] == "EXPOSED", sum(aesifup_input[[eventCol]])]
  n_out_con <- aesifup_input[aesifup_input[[groupCol]] == "CONTROL", sum(aesifup_input[[eventCol]])]
  if (adjust == "unadj") {
    n_pat_exp_weighted <- n_pat_exp
    n_pat_con_weighted <- n_pat_con
    n_out_exp_weighted <- n_out_exp
    n_out_con_weighted <- n_out_con
  }
  else if (adjust == "adj") {
    n_pat_exp_weighted <- aesifup_input[aesifup_input[[groupCol]] == "EXPOSED", sum(aesifup_input[[iptw]])]
    n_pat_con_weighted <- aesifup_input[aesifup_input[[groupCol]] == "CONTROL", sum(aesifup_input[[iptw]])]
    n_out_exp_weighted <- aesifup_input[aesifup_input[[groupCol]] == "EXPOSED", sum(aesifup_input[[eventCol]] * aesifup_input[[iptw]])]
    n_out_con_weighted <- aesifup_input[aesifup_input[[groupCol]] == "CONTROL", sum(aesifup_input[[eventCol]] * aesifup_input[[iptw]])]
  }
  pp_pax <- n_out_exp_weighted/n_pat_exp_weighted * scale_IR
  pp_comp <- n_out_con_weighted/n_pat_con_weighted * scale_IR
  pr <- pp_pax / pp_comp
  pd <- pp_pax - pp_comp

  # estimators - columns to use from boostrap results
  if (!estimate_survival) {
    boot_cols <- c("pp_pax", "pp_comp", "pr", "pd")
    survival_output_list <- list()
  } else {
    boot_cols <- c("rr", "rd", "hr")
    # survival risk estimate: rr, rd, hr
    risk_table <- create_risk_table(
      aesifup = aesifup_input,
      timepoints = risk_window,
      fupCol = "fup",
      eventCol = eventCol,
      use_weights = adjust == "adj",
      iptw = iptw,
      scale_IR = scale_IR,
      comparison_measures = TRUE,
      dummy_code = NA)
    hr_table <- estimate_HR(
      aesifup_input,
      fupCol = "fup",
      eventCol = eventCol,
      iptw = iptw,
      model_type = "crude",
      return_dummy_output = FALSE,
      aesi_name = "",
      dummy_code = NA)
    survival_output_list <- list(
      rr = risk_table$cuminc_est_exp / risk_table$cuminc_est_con,
      rd = risk_table$cuminc_est_exp - risk_table$cuminc_est_con,
      hr = hr_table$hr_est
    )
  }

  # bootstrap CI estimate
  if (is.null(bootstrap))
    return(data.table(n_pat_exp, n_pat_con, n_out_exp, n_out_con,
                      pp_pax, pp_comp, pr, pd, do.call(cbind, survival_output_list)))
  else{
    boot_ci_list <- list()
    for (prev_suffix in boot_cols) {
      boot_col_name <- paste0(prev_suffix, "_", adjust)
      boot_vec <- bootstrap[, get(boot_col_name)]
      boot_ci_list[[prev_suffix]] <- est_inc_prev(
        scale_IR = 1,
        type = "prevalence",
        CImethod = "bootstrap",
        boot_vec = boot_vec)[
          , .(ir_lb, ir_ub)]
    }
    return(data.table(n_pat_exp, n_pat_con, n_out_exp, n_out_con,
                      pp_pax, pp_comp, pr, pd,
                      do.call(cbind, survival_output_list),
                      do.call(cbind, boot_ci_list)))
  }
}
