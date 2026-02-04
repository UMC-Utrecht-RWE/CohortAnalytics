# wrapper function to create table of risks using a KM model
# splits dataset into exposed and control parts
# fits KM model using survfit, grabbing appropriate pyr and eventCount columns
# feeds this model to est_km function for vector of timepoints
# returns 1-KM scaled with scale_IR, with lower and upper bounds
# collects aesi name, time, exposed and control risks into one table
#' @param aesifup data
#' @param timepoints a vector of time to be used in survival model
#' @param fupcol column name of follow up time, to be used in survival model
#' @param eventCol event column name
#' @param use_weights TRUE/FALSE
#' @param iptw column name of inverse probability weight
#' @param scale_IR report rate to multiply by this scale
#' @param comparison_measures TRUE means risk rate will be compared with control group
#' @param dummy_code use code to fill NA values
#' @export

create_risk_table <- function(aesifup,
                              timepoints,
                              fupCol = "fup",
                              eventCol = "eventCount",
                              use_weights = TRUE,
                              iptw = "ip_weight",
                              scale_IR = 10000,
                              comparison_measures = TRUE,
                              dummy_code = -99){
  # aesifup <- aesifup_input_tmp
  # timepoints = max_fuptime
  # target_aesi = target_aesi
  # use_weights = FALSE
  # fupCol = fupCol
  # eventCol = eventCol
  # iptw = iptw
  # scale_IR = scale_IR
  # comparison_measures = comparison_measures

  # if no data, return NA flags
  if(nrow(aesifup) == 0){
    return(data.frame(time =dummy_code,
                      "cuminc_est_exp" = dummy_code,
                      "cuminc_lb_exp"  = dummy_code,
                      "cuminc_ub_exp" = dummy_code,
                      "cuminc_est_con" = dummy_code,
                      "cuminc_lb_con" = dummy_code,
                      "cuminc_ub_con" = dummy_code))
  } else {
    if(comparison_measures){

      # split file into control and exposed part
      dfexp <- copy(aesifup[aesifup$group == "EXPOSED",])
      dfcon <- copy(aesifup[aesifup$group == "CONTROL",])

      if(nrow(dfexp) != 0){
        if(use_weights == TRUE){
          # fit survival model adjusting using weights
          kmexp <- survival::survfit(survival::Surv(get(fupCol), get(eventCol))~1,
                                     cluster = person_id, robust = T, data = dfexp,
                                     weights = get(iptw))
        }else{
          kmexp <- survival::survfit(survival::Surv(get(fupCol), get(eventCol))~1,
                                     cluster = person_id, robust = T, data = dfexp)
        }

        # obtain risks at each timepoint of interest, scaled appropriately
        risk_exp <- as.data.frame(do.call("rbind",
                                          sapply(timepoints, function(s) est_km(kmexp, s, per_pyr = scale_IR),
                                                 simplify = FALSE)))
      } else {
        risk_exp <- data.frame(
          time = dummy_code,
          cuminc_est = 0,
          cuminc_lb = 0,
          cuminc_ub = 0
        )
      }

      # repeat for controls
      if(nrow(dfcon) != 0){
        if(use_weights == TRUE){
          kmcon <-  survival::survfit(survival::Surv(get(fupCol), get(eventCol))~1,
                                      id = person_id, robust = T, data = dfcon,
                                      weights = get(iptw))
        }else{
          kmcon <-  survival::survfit(survival::Surv(get(fupCol), get(eventCol))~1,
                                      id = person_id, robust = T, data = dfcon)
        }
        risk_con <- as.data.frame(do.call("rbind",
                                          sapply(timepoints, function(s) est_km(kmcon, s, per_pyr = scale_IR),
                                                 simplify = FALSE)))
      } else {
        risk_con <- data.frame(
          time = dummy_code,
          cuminc_est = 0,
          cuminc_lb = 0,
          cuminc_ub = 0
        )
      }


    } else if (!comparison_measures){
      # split file into control and exposed part
      dfexp <- copy(aesifup)

      if(use_weights == TRUE){
        # fit survival model adjusting using weights
        kmexp <- survival::survfit(survival::Surv(get(fupCol), get(eventCol))~1,
                                   cluster = person_id, robust = T, data = dfexp,
                                   weights = get(iptw))
      }else{
        kmexp <- survival::survfit(survival::Surv(get(fupCol), get(eventCol))~1,
                                   cluster = person_id, robust = T, data = dfexp)
      }
      # obtain risks at each timepoint of interest, scaled appropriately
      risk_exp <- as.data.frame(do.call("rbind",
                                        sapply(timepoints, function(s) est_km(kmexp, s, per_pyr = scale_IR),
                                               simplify = FALSE)))
    }

    # create output
    output <- data.frame(time = risk_exp$time,
                         "cuminc_est_exp" = risk_exp$cuminc_est,
                         "cuminc_lb_exp"  = risk_exp$cuminc_lb,
                         "cuminc_ub_exp" = risk_exp$cuminc_ub,
                         "cuminc_est_con"  = ifelse(comparison_measures, risk_con$cuminc_est, dummy_code),
                         "cuminc_lb_con" = ifelse(comparison_measures, risk_con$cuminc_lb, dummy_code),
                         "cuminc_ub_con" = ifelse(comparison_measures, risk_con$cuminc_ub, dummy_code))
  }
  return(output)
}


# internal helper function for wrangling km model fits and computing survivals
# function which takes a fitted model from survfit
# computes 1-km for a vector of time points (given in days by default)
# set days_to_years TRUE if model fit objects time is given in years

# if ndaystotimepoint = NULL we just take the last day of follow-up

# outputs time (in days), and point estimates + CIs scaled according to per_pyr
# assumes fit object person time counted in years
est_km <- function(fit, ndaystotimepoint, days_to_years = FALSE, per_pyr = 1){
  #ndaystotimepoint <- 180
  time_out <- ndaystotimepoint
  if(days_to_years == TRUE){
    ndaystotimepoint <- ndaystotimepoint/365.25
  }

  if(!is.null(ndaystotimepoint)){

    if(last(fit$time) >= ndaystotimepoint & min(fit$time) < ndaystotimepoint) { #if not met, we do not show the results
      here <- which(fit$time == max(fit$time[fit$time <= ndaystotimepoint]))
      rsk <- (1 - fit$surv[!is.na(fit$surv)][here])
      rsk.lb <- (1 - fit$upper[!is.na(fit$surv)][here])
      rsk.ub <- (1 - fit$lower[!is.na(fit$surv)][here])
    } else if (min(fit$time) > ndaystotimepoint){
      rsk <- 0
      rsk.lb <- 0
      rsk.ub <- 0
    }
    else {
      rsk <- NA
      rsk.lb <- NA
      rsk.ub <- NA
    }

  } else if (is.null(ndaystotimepoint)){
    time_out <- 9999
    rsk <- (1 - last(fit$surv[!is.na(fit$surv)]))
    rsk.lb <- (1 - last(fit$upper[!is.na(fit$surv)]))
    rsk.ub <- (1 - last(fit$lower[!is.na(fit$surv)]))
  }

  if (rsk.lb < 0 & !is.na(rsk.lb)) rsk.lb <- 0
  if (is.na(rsk.lb)) rsk.lb <- NA


  return(data.frame(time = time_out,
                    cuminc_est =rsk*per_pyr,cuminc_lb = rsk.lb*per_pyr, cuminc_ub = rsk.ub*per_pyr))
}
