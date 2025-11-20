# ---------- dummy output for overall compute cohort -----
# return a single row with optional flag
dummy_compute_rates <- function(target_aesi, end_risk, flag = -99, comparison_measures = TRUE){
  dummy_out <- data.frame(
    "aesi" = target_aesi,
    "n_pat_exp" = flag,
    "n_out_exp" = flag,
    "py_exp" = flag,
    "ir_est_exp" = flag,
    "ir_lb_exp" = flag,
    "ir_ub_exp" = flag,
    "n_pat_con" = flag,
    "n_out_con" = flag,
    "py_con" = flag,
    "ir_est_con" = flag,
    "ir_lb_con" = flag,
    "ir_ub_con" = flag,
    "end_risk" = end_risk,
    "time" = flag,
    "cuminc_est_exp_crude" = flag,
    "cuminc_lb_exp_crude" = flag,
    "cuminc_ub_exp_crude" = flag,
    "cuminc_est_con_crude" = flag,
    "cuminc_lb_con_crude" = flag,
    "cuminc_ub_con_crude" = flag,
    "cuminc_est_exp_adj" = flag,
    "cuminc_lb_exp_adj" = flag,
    "cuminc_ub_exp_adj" = flag,
    "cuminc_est_con_adj" = flag,
    "cuminc_lb_con_adj" = flag,
    "cuminc_ub_con_adj" = flag
  )

  if(comparison_measures){
    dummy_out <- cbind(
      dummy_out,
      data.frame(
        "rr_est_crude" = flag,
        "rr_lb_crude" = flag,
        "rr_ub_crude" = flag,
        "rr_est_adj" = flag,
        "rr_lb_adj" = flag,
        "rr_ub_adj" = flag,
        "rd_est_crude" = flag,
        "rd_lb_crude" = flag,
        "rd_ub_crude" = flag,
        "rd_est_adj" = flag,
        "rd_lb_adj" = flag,
        "rd_ub_adj" = flag,
        "rr_type" = flag,
        "rd_type" = flag
      )
    )
  }

  return(dummy_out)
}
