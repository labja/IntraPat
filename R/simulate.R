#' Simulates multiple trials
#'
#' @param method Vector of dose-finding designs to be employed. The designs will be applied sequentially with the estimated MTD from the previous design as the starting dose. Possible options are "Intra" for Intra-patient escalation, "TPT" for the 3+3 design, "BCRM" for (Bayesian) Continual Reassessment Method and "BOIN" for Bayesian Optimal Interval Design.
#' @param tox_rates Vector of true toxicity rates for each dose level in ascending order
#' @param target Target toxicity rate
#' @param cohortsize Number of patients in a cohort
#' @param start Dose level to be used at the beginning of the trial
#' @param bcrm_constrain Compare bcrm::bcrm
#' @param bcrm_ff Compare bcrm::bcrm
#' @param bcrm_method Compare bcrm::bcrm
#' @param bcrm_prior Compare bcrm::bcrm
#' @param bcrm_stop Compare bcrm::bcrm
#' @param boin_ncohort Compare BOIN::get.oc
#' @param intra_days How many single treatments are given
#' @param intra_split_tox_rates Boolean value whether toxicity rates are split evenly assuming identical and independent toxicity rates for each treatment. FALSE implies a binary toxicity rate.
#' @param tpt_dose_expansion Boolean value whether a dose expansion stage is included if fewer than 6 patients were treated at the estimated MTD
#' @param nsim Number of simulated trials
#' @param seed Seed for simulation
#' @return res Data frame with number of patients, number of DLTs and estimated MTD for each trial
#' @return res_dose Data frame with number of patients and number of DLTs for each combination of dose level and trial
#' @return true_mtd Dose level which is the closest to target toxicity
#'
#' @examples
#'
#' sim_trials(method=c("Intra","TPT"),tox_rates=seq(0.1,0.5,0.1),target=0.3,intra_days=4,nsim=10)

sim_trials <- function(method,tox_rates,target,cohortsize=1,start=1,
                       bcrm_constrain=NULL,bcrm_ff="power",bcrm_method="exact",
                       bcrm_prior=NULL,bcrm_stop=NULL,boin_ncohort=NULL,
                       intra_days=NULL,intra_split_tox_rates=TRUE,
                       tpt_dose_expansion=TRUE,nsim,seed=1) {

  set.seed(seed)

  # Initialize result data frame
  res <- data.frame(n_pat=numeric(nsim),n_dlt=numeric(nsim),mtd_est=numeric(nsim))
  res_dose <- data.frame(dose_level=numeric(),n=numeric(),dlt=numeric())

  true_mtd <- which.min(abs(target-tox_rates))

  # Simulate individual trials and save key numbers
  for (i in 1:nsim) {

    print(paste("Trial",i,sep=" "))

    summary <- summary_trial(sim_trial(method=method,tox_rates=tox_rates,target=target,
                                       cohortsize=cohortsize,start=start,
                                       bcrm_constrain=bcrm_constrain,bcrm_ff=bcrm_ff,bcrm_method=bcrm_method,
                                       bcrm_prior=bcrm_prior,bcrm_stop=bcrm_stop,
                                       boin_ncohort=boin_ncohort,
                                       intra_days=intra_days,intra_split_tox_rates=intra_split_tox_rates,
                                       tpt_dose_expansion=tpt_dose_expansion))

    res$n_pat[i] <- summary$n_pat
    res$n_dlt[i] <- summary$n_dlt
    res$mtd_est[i] <-summary$mtd_est
    res_dose <- rbind(res_dose,dplyr::mutate(summary$res_dose,trial=i))

  }


  res$mtd_est <- tidyr::replace_na(res$mtd_est,0)

  input <- list(method=method,tox_rates=tox_rates,target=target,cohortsize=cohortsize,start=start,
                bcrm_constrain=bcrm_constrain,bcrm_ff=bcrm_ff,bcrm_method=bcrm_method,
                bcrm_prior=bcrm_prior,bcrm_stop=bcrm_stop,boin_ncohort=boin_ncohort,
                intra_days=intra_days,intra_split_tox_rates=intra_split_tox_rates,
                tpt_dose_expansion=tpt_dose_expansion,nsim=nsim,seed=seed)

  return(list(res=res,res_dose=res_dose,true_mtd=true_mtd,input=input))
}


#' Simulates a single trial
#'
#' @param method Vector of dose-finding designs to be employed. The designs will be applied sequentially with the estimated MTD from the previous design as the starting dose. Possible options are "Intra" for Intra-patient escalation, "TPT" for the 3+3 design, "BCRM" for (Bayesian) Continual Reassessment Method and "BOIN" for Bayesian Optimal Interval Design.
#' @param tox_rates Vector of true toxicity rates for each dose level in ascending order
#' @param target Target toxicity rate
#' @param cohortsize Number of patients in a cohort
#' @param start Dose level to be used at the beginning of the trial
#' @param bcrm_constrain Compare bcrm::bcrm
#' @param bcrm_ff Compare bcrm::bcrm
#' @param bcrm_method Compare bcrm::bcrm
#' @param bcrm_prior Compare bcrm::bcrm
#' @param bcrm_stop Compare bcrm::bcrm
#' @param boin_ncohort Compare BOIN::get.oc
#' @param intra_days How many single treatments are given
#' @param intra_split_tox_rates Boolean value whether toxicity rates are split evenly assuming identical and independent toxicity rates for each treatment. FALSE implies a binary toxicity rate.
#' @param tpt_dose_expansion Boolean value whether a dose expansion stage is included if fewer than 6 patients were treated at the estimated MTD
#' @return res Data frame with number of patients, number of DLTs and estimated MTD for each trial
#' @return res_dose Data frame with number of patients and number of DLTs for each combination of dose level and trial
#' @return true_mtd Dose level which is the closest to target toxicity
#'
#' @examples
#'
#' sim_trial(method=c("Intra","TPT"),tox_rates=seq(0.1,0.5,0.1),target=0.3,intra_days=4)

sim_trial <- function(method,tox_rates,target,cohortsize=1,start=1,
                      bcrm_constrain=NULL,bcrm_ff="power",bcrm_method="exact",
                      bcrm_prior=NULL,bcrm_stop=NULL,boin_ncohort=NULL,
                      intra_days=NULL,intra_split_tox_rates=TRUE,
                      tpt_dose_expansion=TRUE) {



  mtd_est <- start
  true_mtd <- which.min(abs(target-tox_rates))
  res_pat <- data.frame(pat=numeric(),dose_level=numeric(),dlt=numeric())

  # Simulate dose-finding trial. Each iteration corresponds to a dose-finding design specified in method.

  for (i in 1:length(method)) {


      res_stage <-  switch(method[i]
                     ,Intra=intra(tox_rates=tox_rates,start=mtd_est,days=intra_days,split_tox_rates=intra_split_tox_rates)
                     ,TPT=tpt(tox_rates=tox_rates,start=mtd_est,dose_expansion=tpt_dose_expansion)
                     ,BCRM=bcrm_jl(tox_rates=tox_rates,target=target,cohortsize=cohortsize,res_pat=res_pat,start=mtd_est,
                                   constrain=bcrm_constrain,ff=bcrm_ff,bcrm_method=bcrm_method,prior=bcrm_prior,bcrm_stop=bcrm_stop)
                     ,BOIN=boin(p.true=tox_rates,target=target,cohortsize=cohortsize,startdose=mtd_est,ncohort=boin_ncohort))
      res_pat <- rbind(res_pat,res_stage$res_pat)
      mtd_est <- res_stage$mtd_est
    }

  input <-  list(method=method,tox_rates=tox_rates,target=target,cohortsize=cohortsize,start=start,
                 bcrm_constrain=bcrm_constrain,bcrm_ff=bcrm_ff,bcrm_method=bcrm_method,
                 bcrm_prior=bcrm_prior,bcrm_stop=bcrm_stop,boin_ncohort=boin_ncohort,
                 intra_days=intra_days,intra_split_tox_rates=intra_split_tox_rates,
                 tpt_dose_expansion=tpt_dose_expansion)

  return(list(res_pat=res_pat,mtd_est=mtd_est,input=input))


}




# Simulates whether DLT occurs or not

# tox_rates Vector of true toxicity rates for each dose level in ascending order
# dose_level Integer corresponding to the element of tox_rates vector
# cohortsize Number of patients in a cohort



sim_dlt <- function(tox_rates,dose_level,cohortsize=1) {


  if (dose_level>length(tox_rates)) stop("Dose Level too high")
  if (dose_level<1) stop("Dose Level too low")

  dlt <- rbinom(n = cohortsize, size = 1, prob = tox_rates[dose_level])
  return(dlt)
}


#' Summarizes the result of a single simulated trial
#'
#' @param res Output from \code{sim_trial}

summary_trial <- function(res) {


  res_pat <- res$res_pat
  mtd_est <- res$mtd_est

  # Get number of patients
  help_df <- res_pat %>%
    dplyr::group_by(design) %>%
    dplyr::slice(which.max(pat))
  n_pat <- sum(help_df$pat)

  # Get number of DLTS
  n_dlt <- sum(res_pat$dlt)


  # Get number of patients by dose level if design is not Intra
  res_dose <- res_pat %>%
    dplyr::filter(design!="Intra") %>%
    dplyr::group_by(dose_level) %>%
    dplyr::summarise(n = dplyr::n(), dlt = sum(dlt), .groups="drop_last")

  # Get number of patients by dose level if design is Intra. Each treatment is counted as 1/intra_days.
    if (any(res_pat$design=="Intra")) {

  days <- res$input$intra_days

  res_dose1 <- res_pat %>%
    dplyr::filter(design=="Intra") %>%
    dplyr::group_by(dose_level) %>%
    dplyr::summarise(n = dplyr::n()/days, dlt = sum(dlt), .groups="drop_last")
    res_dose <- rbind(res_dose1,res_dose)
  }

  # Aggregate dose level for different designs
  res_dose <- res_dose %>%
    dplyr::group_by(dose_level) %>%
    dplyr::summarise_all(list(sum))

  return(list(n_pat=n_pat,n_dlt=n_dlt,mtd_est=mtd_est,res_dose=res_dose))
}

#' Summarizes the result of a multiple simulated trials
#'
#' @param res_list Output from \code{sim_trials}


summary_trials <- function(res_list) {

  true_mtd <- res_list$true_mtd

  # get median number of patients, median number of DLT and median estimated MTD
  res_median <- t(data.frame(apply(res_list$res,2,median)))


  # Get frequency of treatments at dose level
  freq_dose <- res_list$res_dose %>%
    dplyr::group_by(dose_level) %>%
    dplyr::summarise(n_pat=sum(n),.groups="drop") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(freq_pat=n_pat/sum(n_pat),true_mtd=dose_level==true_mtd)

  # Get frequency of RP2D at dose level
  freq_mtd_est <- res_list$res %>%
    dplyr::group_by(mtd_est) %>%
    dplyr::summarise(n_rp2d=dplyr::n(),.groups="drop") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(freq_rp2d=n_rp2d/sum(n_rp2d),true_mtd=mtd_est==true_mtd)

  # Remove this later because included in freq_mtd_est as well?
  accuracy <- sum(res_list$res$mtd_est==res_list$true_mtd)/length(res_list$res$mtd_est)


  return(list(res_median=res_median,accuracy=accuracy,freq_dose=freq_dose,freq_mtd_est=freq_mtd_est))


}
