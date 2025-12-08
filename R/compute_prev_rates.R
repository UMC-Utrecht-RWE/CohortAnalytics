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
  if (adjust == "unadj"){
    #===== SUM PER GROUP
    n_pat_exp <- aesifup_input[group == "EXPOSED", .N]
    n_pat_con <- aesifup_input[group == "CONTROL", .N]
    # sum outcomes
    n_out_exp <- aesifup_input[group == "EXPOSED", sum(.SD[[eventCol]])]
    n_out_con <- aesifup_input[group == "CONTROL", sum(.SD[[eventCol]])]
  } else if (adjust == "adj") {
    n_pat_exp <- aesifup_input[group == "EXPOSED", sum(wt)]
    n_pat_con <- aesifup_input[group == "CONTROL", sum(wt)]
    n_out_exp <- aesifup_input[aesifup_input$group == "EXPOSED", sum(.SD[[eventCol]] * wt)]
    n_out_con <- aesifup_input[aesifup_input$group == "CONTROL", sum(.SD[[eventCol]] * wt)]
  }
  # point estimates
  pp_pax <- n_out_exp / n_pat_exp
  pp_comp <- n_out_con / n_pat_con
  pr <- pp_pax / pp_comp
  pd <- pp_pax - pp_comp

  # CI of pp (pax, unexp) or pr or pd + unweighted or weighted
  boot_ci_list <- list()
  for (prev_suffix in c("pp_pax", "pp_comp", "pr", "pd")) {
    boot_col_name <- paste0(prev_suffix, "_", adjust)
    boot_vec <- bootstrap[, get(boot_col_name)]
    boot_ci_list[[prev_suffix]] <- est_inc_prev(scale_IR = 1, type =  "prevalence", CImethod = "bootstrap", boot_vec = boot_vec)[, .(ir_lb, ir_ub)]
  }
  return(data.table(n_pat_exp, n_pat_con, n_out_exp, n_out_con, pp_pax, pp_comp, pr, pd, do.call(cbind, boot_ci_list)))
}
