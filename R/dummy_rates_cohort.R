# function to create dummy output for the cohort analyses
# to be called if variable missing or no data
# here we must mimic the structure of the output for compute_rates_cohort
dummy_rates_cohort <- function(target_aesi){

  ir_list_exp <- ir_list_con <-
    data.frame(
      ir = -99,
      lb_ir = -99,
      ub_ir = -99)
  # rename incidence rates to be group labelled
  colnames(ir_list_exp) <- paste0(colnames(ir_list_exp),"_exp")
  colnames(ir_list_con) <- paste0(colnames(ir_list_con),"_con")

  ests_ir <- data.frame(aesi = target_aesi,
                        n_pat_exp = -99,
                        n_out_exp = -99,
                        py_exp = -99,
                        ir_list_exp,
                        n_pat_con = -99,
                        n_out_con = -99,
                        py_con = -99,
                        ir_list_con)

  ests_crude <- ests_adj <- data.frame(irr_est = -99,
                                       lb_irr = -99,
                                       ub_irr = -99,
                                       ird = -99,
                                       lb_ird = -99,
                                       ub_ird = -99)

  # add target_aesi to crude and adjusted estimates
  ests_crude$aesi <- ests_adj$aesi <- target_aesi

  # return a list with three objects to be combined later
  return(list(ests_ir = ests_ir,
              ests_crude = ests_crude,
              ests_adj = ests_adj))
}
