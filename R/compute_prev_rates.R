#' Wrapper function to compute statistics relating to prevelance
#'
#'  rates, ratios and differences
#'  per group and between groups
#'  with or without weights
#'
#' key functionality:
#'  est_inc_prev
#'
#' output:
#'  point estimates (mean)
#'  CI based on bootstrap results
#'
#' @param aesifup_input data
#' @param eventCol column of event 0/1
#' @param iptw column of weight
#' @param target_aesi choice of outcome to calculate prevlance
#' @param bootstrap bootstrap result to calculate CI
#' @param adjust either 'adj' adjusted or 'unadj' unadjusted by iptw
#'
#' @export
compute_prev_rates <- function(
    aesifup_input,
    eventCol = "event",
    iptw = "wt",
    target_aesi,
    bootstrap,
    adjust = "adj"
) {
  n_pat_exp <- aesifup_input[group == "EXPOSED", .N]
  n_pat_con <- aesifup_input[group == "CONTROL", .N]
  n_out_exp <- aesifup_input[group == "EXPOSED", sum(.SD[[eventCol]])]
  n_out_con <- aesifup_input[group == "CONTROL", sum(.SD[[eventCol]])]
  if (adjust == "unadj") {
    n_pat_exp_weighted <- n_pat_exp
    n_pat_con_weighted <- n_pat_con
    n_out_exp_weighted <- n_out_exp
    n_out_con_weighted <- n_out_con
  }
  else if (adjust == "adj") {
    n_pat_exp_weighted <- aesifup_input[group == "EXPOSED", sum(wt)]
    n_pat_con_weighted <- aesifup_input[group == "CONTROL", sum(wt)]
    n_out_exp_weighted <- aesifup_input[group == "EXPOSED", sum(.SD[[eventCol]] *
                                                                  wt)]
    n_out_con_weighted <- aesifup_input[group == "CONTROL", sum(.SD[[eventCol]] *
                                                                  wt)]
  }
  pp_pax <- n_out_exp_weighted/n_pat_exp_weighted * scale_IR
  pp_comp <- n_out_con_weighted/n_pat_con_weighted * scale_IR
  pr <- pp_pax / pp_comp
  pd <- pp_pax - pp_comp
  if(is.null(bootstrap))
    return(data.table(n_pat_exp, n_pat_con, n_out_exp, n_out_con,
                      pp_pax, pp_comp, pr, pd))
  else{
    boot_ci_list <- list()
    for (prev_suffix in c("pp_pax", "pp_comp", "pr", "pd")) {
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
                      pp_pax, pp_comp, pr, pd, do.call(cbind, boot_ci_list)))
  }
}
