# wrapper function to create outputs of incidence rates, IRRs, IRDs adjusted and crude
# assumes a single aesi cohort input file, with columns eventCount_<analysis_id> and pyr_<analysis_id>
# Computes and outputs (for eeach control and exposed group)
# variablename
# system
# aesi label
# compute; n_patients, n_outcomes, person-years
# ir per 100,000 py, 95 CI
# unless option is prevalence, in which case ....
# Crude IRR, 95 CI
# Crude IRD, 95 CI
# Adjusted IRR, 95 CI
# Adjusted IRD, 95 CI

#' Wrapper function to compute statistics relating to incidence, prevalences and risks (ratios and differences)
#' per group and between groups
#'
#' @param aesifup_input data.table, typically representing two groups (exposed and control), with columns specifying person-time, occurence of events, and group membership. Must contain a column `group` with values EXPOSED and CONTROL. May be outputted by `process_cohort_data`
#' @param fupCol name of column representing days of follow-up time. defaults to "fup"
#' @param pyrCol name of column representing person-years of follow-up time. defaults to "pyr"
#' @param pyr_offset name of column which is used as input to poisson (GEE) models if necessary, typically person-time. "pyr_offset"
#' @param eventCol name column identifying binary event count for this unit, "eventCount"
#' @param iptw name column containing (ipt or other) weights to be used in the analysis
#' @param target_aesi character string giving name of the outcome variable, for informative logging
#' @param type "incidence" or "prevalence", defining if incidence rates or prevalenes are to be computed
#' @param CImethod defaults to "wilson"; option passed to `est_inc_prev` defining formula for analytic CI computation
#' @param risk_type defaults "survival". if set to survival, function outputs 1-KM risk differences and ratios. Otherwise prevalence/incidence ratios and differences, as specified by `type`
#' @param CumulativeInc TRUE/FALSE, compute cumulative incidences using 1-KM. If not, skip and return -99 flag
#' @param comparison_measures controls whether 1-KM is based on the latest available follow-up in both exposed and control (TRUE) or just the maximum follow-up in the dataset as a whole (FALSE)
#' @param scale_IR what scaling to perform for incidence/prevalence, i.e., events per X. Defaults to 10000
#' @param model_based_control TRUE/FALSE. Use GEE to estimate incidence/prevalence in the control group (TRUE) or an analytic formula that assumes independent observations (FALSE)
#' @param weighted_IR use iptw column to weight the incidence/prevalence estimates per group (TRUE) or output unweighted incidence/prevalence in the exposed and control (FALSE)
#' @param end_risk
#' @param output_format defaults to `data.table`, otherwise returns data.frame.
#' @param minimum_count_for_comparative defaults to 3; should a minimum event count be applied in order to display results, all estimates relating to event counts less than this will be suppressed
#'
#' @returns
#' @export
#'
#' @examples
compute_rates_cohort <- function(aesifup_input,
                                 fupCol = "fup",
                                 pyrCol = "pyr",
                                 pyr_offset = "pyr_offset",
                                 eventCol = "eventCount",
                                 iptw = "ip_weight",
                                 target_aesi,
                                 type = "incidence",
                                 CImethod = "wilson",
                                 risk_type = "survival",
                                 CumulativeInc = TRUE,
                                 comparison_measures = TRUE,
                                 scale_IR = 10000,
                                 model_based_control = FALSE,
                                 weighted_IR = FALSE,
                                 end_risk,
                                 output_format = "data.table",
                                 minimum_count_for_comparative = 3){

  # aesifup_input = aesifup_input_tmp
  # # global settings
  # fupCol
  # pyrCol
  # pyr_offset
  # eventCol
  # iptw
  # model_based_control = FALSE
  # comparison_measures = comparison_measures
  # CumulativeInc = CumulativeInc
  # # aesi / analysis-specific settings
  # target_aesi = target_aesi
  # type = estimator #incidence/prevalence
  # risk_type = risk_type # based on km or incidence/prevalence
  # end_risk = end_risk


  # if no data, return NA flags
  # if(nrow(aesifup_input) == 0){
  #   # return(data.frame(VariableName = target_aesi, NegativeControlOutcome = label, time = risk_exp$time,
  #   #                   "km_rsk_exp" =  -88, "lb_km_exp" = -88, "ub_km_exp" = -88,
  #   #                   "km_rsk_con" =  -88, "lb_km_con" = -88, "ub_km_con" = -88))
  # }
  #====================================   COUNTS    =======================================#
  #===== SUM PER GROUP
  n_pat_exp <- aesifup_input[aesifup_input$group == "EXPOSED", .N]
  n_pat_con <- aesifup_input[aesifup_input$group == "CONTROL", .N]

  # sum outcomes
  n_out_exp <- aesifup_input[aesifup_input$group == "EXPOSED", sum(.SD[[eventCol]])]
  n_out_con <- aesifup_input[aesifup_input$group == "CONTROL", sum(.SD[[eventCol]])]

  # person-years
  py_exp <- aesifup_input[aesifup_input$group == "EXPOSED", sum(.SD[[pyrCol]])]
  py_con <- aesifup_input[aesifup_input$group == "CONTROL", sum(.SD[[pyrCol]])]

  #==================================  INCIDENCE RATES   ==================================#
  # create appropriate model type depending on desired estimand (incidence/prevalence)
  # with GEE model-based approahces incidence vs prevalence is determined by the offset
  # incidence - offset = log(pyr)
  # prevalence - offset = 1
  model_formula <- as.formula(paste0(eventCol, " ~ group + offset(",pyr_offset,")"))
  model_formula_no_group <- as.formula(paste0(eventCol, " ~ offset(",pyr_offset,")"))

  # Fit crude and adjusted models to extract incidence/prevalence rate in each group
  # relies on internal function defined below in this script
  # if valid fit, returns fit object
  # if not, returns character string, prints to log
  # aesifup must have person_id_num, iptw for adjusted model
  model_crude <- fitmod_gee(model_formula, model_type = "crude", aesi_name = target_aesi,
                            aesifup_input = aesifup_input)
  model_adj <- fitmod_gee(model_formula, model_type = "adj", iptw = iptw, aesi_name = target_aesi,
                          aesifup_input = aesifup_input)
  # estimate prevalence ratios and

  # ---- incidence rate / prevalence proportion statistics -----
  # incidence rates /prevalence proportiona and CI
  # computed analytically or with help of a model (clustering)

  # compute rates for EXPOSED
  # if weighted IR requested, always use model based IR for exposed
  if(weighted_IR == TRUE){
    model_adj_exposed <- fitmod_gee(model_formula_no_group,
                                    model_type = "adj",
                                    iptw = iptw,
                                    aesi_name = target_aesi,
                                    aesifup_input = aesifup_input[aesifup_input$group == "EXPOSED", ])
    ir_list_exp <- est_inc_prev_model(model_adj_exposed, group = "EXPOSED", scale_IR = scale_IR)

  } else {
    ir_list_exp <- est_inc_prev(n_pat_exp, n_out_exp, py_exp,scale_IR = scale_IR, type = type, CImethod = CImethod)
  }

  # compute rates for CONTROL
  if(model_based_control){
    if(weighted_IR == TRUE){
      model_adj_control <- fitmod_gee(model_formula_no_group,
                                      model_type = "adj",
                                      iptw = iptw,
                                      aesi_name = target_aesi,
                                      aesifup_input = aesifup_input[aesifup_input$group == "CONTROL", ])
      ir_list_con <- est_inc_prev_model(model_adj_control, group = "CONTROL", scale_IR = scale_IR)
    }else{
      model_crude_control <- fitmod_gee(model_formula_no_group,
                                        model_type = "crude",
                                        aesi_name = target_aesi,
                                        aesifup_input = aesifup_input[aesifup_input$group == "CONTROL", ])
      ir_list_con <- est_inc_prev_model(model_crude_control,group = "CONTROL", scale_IR = scale_IR)
    }
  } else{
    ir_list_con <- est_inc_prev(n_pat_con, n_out_con, py_con,scale_IR = scale_IR, type = type, CImethod = CImethod)
  }

  # rename incidence rates to be group labelled
  colnames(ir_list_exp) <- paste0(colnames(ir_list_exp),"_exp")
  colnames(ir_list_con) <- paste0(colnames(ir_list_con),"_con")



  #===============================  CUMULATIVE INCIDENCE   ==================================#
  if(risk_type == "survival" & CumulativeInc){
    # ------------------- compute cumulative incidences/risks -------------------
    # first compute the max follow up time in each group

    if(comparison_measures){
      max_fuptime_tmp <- suppressWarnings(c(aesifup_input[group == "EXPOSED", max(get(fupCol))],
                                            aesifup_input[group == "CONTROL", max(get(fupCol))]))
      max_fuptime <- min(max_fuptime_tmp[!is.infinite(max_fuptime_tmp)])

    } else if(!comparison_measures){
      max_fuptime <- aesifup_input[, max(get(fupCol))]

    }

    # evaluate risks at the max follow up time

    #  matched and weighted
    cuminc_ests_crude <- create_risk_table(aesifup_input,
                                           timepoints = max_fuptime,
                                           target_aesi = target_aesi,
                                           use_weights = FALSE,
                                           fupCol = fupCol,
                                           eventCol = eventCol,
                                           iptw = iptw,
                                           scale_IR = scale_IR,
                                           comparison_measures = comparison_measures)
    if(!is.null(iptw)){
      cuminc_ests_adj <- create_risk_table(aesifup_input,
                                           timepoints = max_fuptime,
                                           target_aesi = target_aesi,
                                           use_weights = TRUE,
                                           fupCol = fupCol,
                                           eventCol = eventCol,
                                           iptw = iptw,
                                           scale_IR = scale_IR,
                                           comparison_measures = comparison_measures)
    } else {
      cuminc_ests_adj <- data.frame(time =-99,
                                    "cuminc_est_exp" = -99,
                                    "cuminc_lb_exp"  = -99,
                                    "cuminc_ub_exp" = -99,
                                    "cuminc_est_con" = -99,
                                    "cuminc_lb_con" = -99,
                                    "cuminc_ub_con" = -99)

    }
  } else { #prevalence
    cuminc_ests_crude <- data.frame(time =-99,
                                    "cuminc_est_exp" = -99,
                                    "cuminc_lb_exp"  = -99,
                                    "cuminc_ub_exp" = -99,
                                    "cuminc_est_con" = -99,
                                    "cuminc_lb_con" = -99,
                                    "cuminc_ub_con" = -99)
    cuminc_ests_adj <- data.frame(time =-99,
                                  "cuminc_est_exp" = -99,
                                  "cuminc_lb_exp"  = -99,
                                  "cuminc_ub_exp" = -99,
                                  "cuminc_est_con" = -99,
                                  "cuminc_lb_con" = -99,
                                  "cuminc_ub_con" = -99)
  }

  # ------ collect output --------
  # add label to cumulative incidence objects
  colnames(cuminc_ests_crude)[-1] <- paste0(colnames(cuminc_ests_crude)[-1],"_crude")
  colnames(cuminc_ests_adj) <- paste0(colnames(cuminc_ests_adj),"_adj")

  # drop time column duplication
  cuminc_ests <- data.frame(cuminc_ests_crude,cuminc_ests_adj[-1])

  # replace NA with non-estimable flag
  cuminc_ests[is.na(cuminc_ests)] <- -88


  if(!comparison_measures) {
    ests_out <- data.frame(aesi = target_aesi,
                           # descriptive estimates relating to exposed
                           n_pat_exp = n_pat_exp,
                           n_out_exp = n_out_exp,
                           py_exp = py_exp,
                           ir_list_exp,
                           # descriptive estiamtes contorl
                           n_pat_con = n_pat_con,
                           n_out_con = n_out_con,
                           py_con = py_con,
                           ir_list_con,
                           # cumulative incidences
                           end_risk = ifelse(is.null(end_risk),"NULL",end_risk),
                           cuminc_ests
    )

  } else if (comparison_measures) {
    #================================   COMPARISION MEAURES   =================================#
    # ----------------------------------------------------------------------------
    # ----------------------- Risk Difference and Ratios -------------------------
    # ----------------------------------------------------------------------------

    #  ------------- for survival-type outcomes, calculate ---------------
    # Hazard Ratios (matched and adjusted)
    # Risk Differences based on cumulative incidence differences
    if(nrow(aesifup) == 0){
      #ratio
      rr_list <- data.frame(rr_est_crude = -99,
                            rr_lb_crude = -99,
                            rr_ub_crude = -99,
                            rr_est_adj = -99,
                            rr_lb_adj = -99,
                            rr_ub_adj = -99)
      #difference
      rd_list <- data.frame(rd_est_crude = -99,
                            rd_lb_crude = -99,
                            rd_ub_crude = -99,
                            rd_est_adj = -99,
                            rd_lb_adj = -99,
                            rd_ub_adj = -99)

      risk_ratio_diff <- data.frame(rr_list,
                                    rd_list,
                                    rr_type = "hazard",
                                    rd_type = "km")
    } else {

      if(risk_type == "survival"){#--------------------------- SURVIVAL
        #---------ratio
        if(n_out_con == 0 | n_out_exp == 0) {
          rr_list <- data.frame(rr_est_crude = -88,
                                rr_lb_crude = -88,
                                rr_ub_crude = -88,
                                rr_est_adj = -88,
                                rr_lb_adj = -88,
                                rr_ub_adj = -88)
        } else {
          crude_hr  <- estimate_HR(aesifup = aesifup_input, fupCol = fupCol,
                                   eventCol = eventCol,
                                   iptw = iptw,
                                   model_type = "crude",
                                   aesi_name = target_aesi)

          adj_hr  <- estimate_HR(aesifup = aesifup_input, fupCol = fupCol,
                                 eventCol = eventCol,
                                 iptw = iptw,
                                 model_type = "adj",
                                 aesi_name = target_aesi)
          #rename columns
          colnames(adj_hr) <- paste0(colnames(adj_hr),"_adj")
          colnames(crude_hr) <- paste0(colnames(crude_hr),"_crude")

          # bind together
          rr_list <- data.frame(cbind(crude_hr, adj_hr))
          colnames(rr_list) <- gsub("hr", "rr", colnames(rr_list))
        }
        #----------difference
        # risk difference - take difference in risks
        rd_est_crude <- cuminc_ests_crude$cuminc_est_exp - cuminc_ests_crude$cuminc_est_con
        rd_est_adj <- cuminc_ests_adj$cuminc_est_exp  - cuminc_ests_adj$cuminc_est_con

        if(is.na(rd_est_crude)) rd_est_crude <- -88
        if(is.na(rd_est_adj)) rd_est_adj <- -88

        if(!is.null(minimum_count_for_comparative)){
          if(n_out_exp < minimum_count_for_comparative | n_out_con < minimum_count_for_comparative){
            rd_list <- data.frame(rd_est_crude = -77,
                                  rd_lb_crude = -77,
                                  rd_ub_crude = -77,
                                  rd_est_adj = -77,
                                  rd_lb_adj = -77,
                                  rd_ub_adj = -77)
          } else {
            rd_list <- data.frame(rd_est_crude = rd_est_crude,
                                  rd_lb_crude = -88, #To be plugged in when the bootstrapping is ready
                                  rd_ub_crude = -88,
                                  rd_est_adj = rd_est_adj,
                                  rd_lb_adj = -88,
                                  rd_ub_adj = -88)
          }
        } else {
          rd_list <- data.frame(rd_est_crude = rd_est_crude,
                                rd_lb_crude = -88, #To be plugged in when the bootstrapping is ready
                                rd_ub_crude = -88,
                                rd_est_adj = rd_est_adj,
                                rd_lb_adj = -88,
                                rd_ub_adj = -88)
        }
        risk_ratio_diff <- data.frame(rr_list,
                                      rd_list,
                                      rr_type = "hazard",
                                      rd_type = "km")


      } else if(risk_type != "survival"){
        # ------- For non-survival type outcomes (incidence/prevalence) -------------
        # compute incidence rate or prevalence proportion ratio
        # compute incidence rate or prevalence proportion difference

        # use internal function defined below
        # this function either processes the model, or returns dummy output if model input is character string
        risk_crude <- extract_ests_gee(model_crude,scale_IR = scale_IR)
        risk_adj <- extract_ests_gee(modelobj = model_adj, scale_IR = scale_IR)

        # change names of columns for consistency
        colnames(risk_crude) <- paste0(colnames(risk_crude),"_crude")
        colnames(risk_adj) <- paste0(colnames(risk_adj),"_adj")

        # colnames(risk_crude) <- gsub("irr", "rr", colnames(risk_crude))
        # colnames(risk_crude) <- gsub("irr", "rr", colnames(risk_adj))

        if(!is.null(minimum_count_for_comparative)){
          if(n_out_exp < minimum_count_for_comparative | n_out_con < minimum_count_for_comparative){
            risk_crude[,4:6] <- -77
            risk_adj[,4:6] <- -77
          }
        }

        # put columns together in correct order
        risk_ratio_diff  <- data.frame(
          # risk ratios
          risk_crude[,1:3],
          risk_adj[,1:3],
          # risk differences
          risk_crude[,4:6],
          risk_adj[,4:6],
          rr_type = type,
          rd_type = type
        )
        colnames(risk_ratio_diff) <- gsub("irr", "rr", colnames(risk_ratio_diff))
        colnames(risk_ratio_diff) <- gsub("ird", "rd", colnames(risk_ratio_diff))

      }
    } # close if(nrow(aesifup) == 0)

    #===================================   COMBINE ALL   ======================================#
    #  ------ collect base information in a data.frame -----
    ests_out <- data.frame(aesi = target_aesi,
                           # descriptive estimates relating to exposed
                           n_pat_exp = n_pat_exp,
                           n_out_exp = n_out_exp,
                           py_exp = py_exp,
                           ir_list_exp,
                           # descriptive estiamtes contorl
                           n_pat_con = n_pat_con,
                           n_out_con = n_out_con,
                           py_con = py_con,
                           ir_list_con,
                           # cumulative incidences
                           end_risk = ifelse(is.null(end_risk),"NULL",end_risk),
                           cuminc_ests,
                           # risk ratios and differences
                           risk_ratio_diff
    )
  } # close if(comparison_measures)
  if(output_format == "data.table") {
    return(as.data.table(ests_out))
  } else {
    # return a list with three objects to be combined later
    return(ests_out)
  }
} # close the function
