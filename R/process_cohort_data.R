# ---- function to perform standard pre-processing steps
process_cohort_data <- function(aesifup,
                                estimator,
                                fupdays = "fup",
                                pyrCol = "pyr",
                                ipweight = "ip_weight") {
  # check if person years in dataset, if not, compute
  if (!(pyrCol %in% names(aesifup))) {
    aesifup[, (pyrCol) := get(fupdays) / 365.25]
  }
  # 1. set CONTROL group as reference category so model parameters always reflect exposure
  set(aesifup, j = "group", value = as.factor(aesifup[["group"]]))
  set(aesifup,
      j = "group",
      value = relevel(aesifup[["group"]], ref = "CONTROL"))
  # 2. create numeric person id column
  # necessary for proper functioning of geepack::geeglm()
  aesifup[, person_id_num := as.numeric(as.factor(person_id))]
  # 3. sort by numeric person id (necessary for geeM)
  setorder(aesifup, person_id_num)
  # create offset for modeling incidence/prevalence
  if (estimator == "prevalence") {
    # 4a. if using a prevalence estimator, create an offest=1 column in the data.table
    aesifup[, pyr_offset := log(1)]
  } else{
    # 4b. create log person years for incidence estimators
    aesifup[, pyr_offset := log(get(pyrCol))]
  }
  #5. check for individuals with exact 0 follow-up, remove if necessary
  zerofup_indicator <- aesifup[, get(pyrCol)] == 0
  if (sum(zerofup_indicator) != 0) {
    # check what reason they have zero follow up time
    # grab pair ids for removal
    drop_pair_ids <- unique(aesifup[zerofup_indicator, id])
    # drop these people and write to log withinformation
    logger::log_info(
      paste(
        "For AESI",
        target_aesi,
        sum(zerofup_indicator),
        "individuals with zero fup... Dropping",
        length(drop_pair_ids),
        "pairs from the analysis"
      )
    )
    aesifup <- aesifup[!id %in% drop_pair_ids, ]
  }

  # 6. check for NA propensity scores, drop
  zeroweight_indicator <- aesifup[,is.na(get(ipweight))]
  if (sum(zeroweight_indicator) != 0) {
    # check what reason they have zero follow up time
    # grab pair ids for removal
    drop_pair_ids_weight <- unique(aesifup[zeroweight_indicator, id])
    # drop these people and write to log withinformation
    logger::log_info(
      paste(
        "For AESI",
        target_aesi,
        sum(zeroweight_indicator),
        "individuals with zero fup... Dropping",
        length(drop_pair_ids_weight),
        "pairs from the analysis"
      )
    )
    aesifup <- aesifup[!id %in% drop_pair_ids_weight, ]
  }

  return(aesifup)
}
