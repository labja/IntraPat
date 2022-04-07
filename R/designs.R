# Simulates intra-patient escalation stage
#
# tox_rates Vector of true toxicity rates for each dose level in ascending order
# start Starting dose level
# days How many single treatments are given
# split_tox_rates  Boolean value whether toxicity rates are split evenly assuming identical and independent toxicity rates for each treatment. FALSE implies a binary toxicity rate.


intra <- function(tox_rates,start,days,split_tox_rates) {

  if (missing(start)) stop("Starting dose level must be provided")

  # Initialization
  if (split_tox_rates) {
    tox_rates <- 1-((1-tox_rates)^(1/days))
  }

  mtd_est <- start
  res_pat <- data.frame(pat=numeric(),day=numeric(),dose_level=numeric(),dlt=logical())
  pat <- 1
  stop <- FALSE

  while (!stop) {


    for (day in 1:days) {


      # Determining dose level

      if (day==1) {

        # Day 1 dose is either previous mtd_est - 1 or 1
        dose_level <- max(mtd_est,2) - 1

      } else {

        # Dose escalation
        dose_level <- mtd_est + 1

      }


      # Simulating DLT and saving results
      dlt <- sim_dlt(tox_rates=tox_rates,dose_level=dose_level)
      new_res <- data.frame(pat=pat,day=day,dose_level=dose_level,dlt=dlt)
      res_pat <- rbind(res_pat,new_res)

      if (dlt) {

        # DLT -> Stop intra-patient escalation
        stop <- TRUE
        mtd_est <- dose_level - 1
        break

      } else {

        mtd_est <- dose_level

        if (mtd_est==length(tox_rates)) {

          stop <- TRUE
          break
        }

      }
    }

    pat <- pat + 1
  }

  res_pat <- dplyr::select(res_pat,-day)
  res_pat$design <- "Intra"

  # Avoid mtd_est below 1
  mtd_est <- max(mtd_est,1)

  return(list(res_pat=res_pat,mtd_est=mtd_est))
}



# Simulates 3+3 stage

# tox_rates Vector of true toxicity rates for each dose level in ascending order
# start Starting dose level
# dose_expansion Boolean value whether or not a dose expansion stage is included if fewer than 6 patients were treated at the estimated MTD

tpt <- function(tox_rates,start,dose_expansion) {

  if (missing(start)) stop("Starting dose level must be provided")

  # Initialization
  mtd_est <- start
  res_pat <- data.frame(pat=numeric(),dlt=numeric(),dose_level=numeric())
  stop <- FALSE

  while (!stop) {


    # Cohort of three patient
    pat1 <- sim_dlt(tox_rates=tox_rates,dose_level=mtd_est)
    pat2 <- sim_dlt(tox_rates=tox_rates,dose_level=mtd_est)
    pat3 <- sim_dlt(tox_rates=tox_rates,dose_level=mtd_est)
    new_pat <- data.frame(pat=(nrow(res_pat)+1):(nrow(res_pat)+3),dlt=c(pat1,pat2,pat3),dose_level=mtd_est)
    (res_pat <- rbind(res_pat,new_pat))

    # Check if trial has to stop and get next dose
    (stop_decision <- tpt_stopcheck(res_pat=res_pat,mtd_est=mtd_est,tox_rates=tox_rates))
    mtd_est <- stop_decision$mtd_est
    stop <- stop_decision$stop


    # Possibility of dose expansion if there were fewer than 6 pat treated at mtd_est
    if (dose_expansion & stop) {

      prev_pat <- res_pat %>%
        dplyr::filter(dose_level==mtd_est)

      if (nrow(prev_pat) == 3 & mtd_est!= 0) {

        pat4 <- sim_dlt(tox_rates=tox_rates,dose_level=mtd_est)
        pat5 <- sim_dlt(tox_rates=tox_rates,dose_level=mtd_est)
        pat6 <- sim_dlt(tox_rates=tox_rates,dose_level=mtd_est)
        new_pat <- data.frame(pat=(nrow(res_pat)+1):(nrow(res_pat)+3),dose_level=mtd_est,dlt=c(pat4,pat5,pat6))
        res_pat <- rbind(res_pat,new_pat)

        stop_decision <- tpt_stopcheck(res_pat=res_pat,mtd_est=mtd_est,tox_rates=tox_rates)
        mtd_est <- stop_decision$mtd_est
        stop <- stop_decision$stop
    }
    }
  }

  res_pat$design <- "TPT"

  return(list(res_pat=res_pat,mtd_est=mtd_est))
}

# Simulates CRM stage

# tox_rates Vector of true toxicity rates for each dose level in ascending order
# target Target toxicity rate
# cohortsize Number of patients in a cohort
# start Starting dose level
# res_pat Toxixity data from previous stage (if there was a previous stage)
# bcrm_method Compare method in bcrm::bcrm
# constrain Compare bcrm::bcrm
# ff Compare bcrm::bcrm
# prior Compare p.tox0 in bcrm::bcrm
# prior.alpha Compare bcrm::bcrm
# bcrm_stop Compare stop in bcrm::bcrm

bcrm_jl <- function(tox_rates,target,cohortsize,start,res_pat,bcrm_method,constrain,ff,
                    prior,prior.alpha=c(1,1,1),bcrm_stop) {



  if (missing(start)) stop("Starting dose level must be provided")

  mtd_est <- start

  # Use previous toxicity data (if previous stage wasn't intra-patient escalation)
  res_pat <- res_pat[res_pat$design!="Intra",]

  stop <- FALSE
  pat <- 1:cohortsize

  while (!stop) {

    dlt <- sim_dlt(tox_rates=tox_rates,dose_level=mtd_est,cohortsize=cohortsize)
    new_res <- data.frame(pat=pat,dose_level=mtd_est,dlt=dlt,design="BCRM")
    res_pat <- rbind(res_pat,new_res)

    # Create a data frame which is compatible with bcrm package from res_pat
    data <- data.frame(patient=1:nrow(res_pat),dose=res_pat$dose_level,tox=res_pat$dlt)


    decision <- bcrm_ndose(stop=bcrm_stop,data=data,p.tox0=prior,target.tox=target,
                       prior.alpha=prior.alpha,ff=ff,constrain=constrain,
                       method=bcrm_method)

    mtd_est <- decision$ndose
    stop <- decision$stop

    pat <- pat + cohortsize

  }

  # If there was a previous stage only return res_pat from BCRM stage
  res_pat <- res_pat[res_pat$design=="BCRM",]

  return(list(res_pat=res_pat,mtd_est=mtd_est))
}


# Simulates BOIN stage

# target Compare BOIN::get.oc
# p.true Compare BOIN::get.oc
# ncohort Compare BOIN::get.oc
# cohortsize Number of patients in a cohort


boin <- function (target, p.true, ncohort, cohortsize, n.earlystop = 100,
                     startdose, titration = FALSE, p.saf = 0.6 * target,
                     p.tox = 1.4 * target, cutoff.eli = 0.95, extrasafe = FALSE,
                     offset = 0.05, ntrial=1) {

  # Most of the code is copied from BOIN package

  if (target < 0.05) {
    stop("the target is too low!")
  }
  if (target > 0.6) {
    stop("the target is too high!")
  }
  if ((target - p.saf) < (0.1 * target)) {
    stop("the probability deemed safe cannot be higher than or too close to the target!")
  }
  if ((p.tox - target) < (0.1 * target)) {
    stop("the probability deemed toxic cannot be lower than or too close to the target!")
  }
  if (offset >= 0.5) {
    stop("the offset is too large!")
  }
  if (n.earlystop <= 6) {
    warning("the value of n.earlystop is too low to ensure good operating characteristics. Recommend n.earlystop = 9 to 18.")
  }

  if (cohortsize == 1)
    titration = FALSE
  lambda_e = log((1 - p.saf)/(1 - target))/log(target * (1 -
                                                           p.saf)/(p.saf * (1 - target)))
  lambda_d = log((1 - target)/(1 - p.tox))/log(p.tox * (1 -
                                                          target)/(target * (1 - p.tox)))
  ndose = length(p.true)
  npts = ncohort * cohortsize
  Y = matrix(rep(0, ndose * ntrial), ncol = ndose)
  N = matrix(rep(0, ndose * ntrial), ncol = ndose)
  dselect = rep(0, ntrial)
  ft = TRUE
  if (cohortsize > 1) {
    temp = BOIN::get.boundary(target, ncohort, cohortsize, n.earlystop = ncohort *
                          cohortsize, p.saf, p.tox, cutoff.eli, extrasafe)$full_boundary_tab
  }
  else {
    temp = BOIN::get.boundary(target, ncohort, cohortsize, n.earlystop = ncohort *
                          cohortsize, p.saf, p.tox, cutoff.eli, extrasafe)$boundary_tab
  }
  b.e = temp[2, ]
  b.d = temp[3, ]
  b.elim = temp[4, ]
  for (trial in 1:ntrial) {
    y <- rep(0, ndose)
    n <- rep(0, ndose)
    earlystop = 0
    d = startdose
    elimi = rep(0, ndose)
    if (titration) {
      z <- (runif(ndose) < p.true)
      if (sum(z) == 0) {
        d = ndose
        n[1:ndose] = 1
      }
      else {
        d = which(z == 1)[1]
        n[1:d] = 1
        y[d] = 1
      }
    }
    for (i in 1:ncohort) {
      if (titration & n[d] < cohortsize && ft) {
        ft = FALSE
        if (d > 1)
          d = d - 1
        y[d] = y[d] + sum(runif(cohortsize - 1) < p.true[d])
        n[d] = n[d] + cohortsize - 1
      }
      else {
        newcohort = runif(cohortsize) < p.true[d]
        if ((sum(n) + cohortsize) >= npts) {
          nremain = npts - sum(n)
          y[d] = y[d] + sum(newcohort[1:nremain])
          n[d] = n[d] + nremain
          break
        }
        else {
          y[d] = y[d] + sum(newcohort)
          n[d] = n[d] + cohortsize
        }
      }
      if (!is.na(b.elim[n[d]])) {
        if (y[d] >= b.elim[n[d]]) {
          elimi[d:ndose] = 1
          if (d == 1) {
            earlystop = 1
            break
          }
        }
        if (extrasafe) {
          if (d == 1 && n[1] >= 3) {
            if (1 - pbeta(target, y[1] + 1, n[1] - y[1] +
                          1) > cutoff.eli - offset) {
              earlystop = 1
              break
            }
          }
        }
      }
      if (n[d] >= n.earlystop && ((y[d] > b.e[n[d]] &&
                                   y[d] < b.d[n[d]]) || (d == 1 && y[d] >= b.d[n[d]]) ||
                                  ((d == ndose || elimi[d + 1] == 1) && y[d] <=
                                   b.e[n[d]])))
        break
      if (y[d] <= b.e[n[d]] && d != ndose) {
        if (elimi[d + 1] == 0)
          d = d + 1
      }
      else if (y[d] >= b.d[n[d]] && d != 1) {
        d = d - 1
      }
      else {
        d = d
      }
    }
    Y[trial, ] = y
    N[trial, ] = n

    # Extract the values needed to return from this function
    dose_level <- unlist(lapply(c(1:length(p.true)),function(x,n){rep(x,n[x])},n=n))
    dlt <- unlist(lapply(c(1:length(p.true)),function(x,n,y){c(rep(0,n[x]-y[x]),rep(1,y[x]))},n=n,y=y))
    pat <- 1:length(dose_level)
    res_pat <- data.frame(pat=pat,dose_level=dose_level,dlt=dlt,design="BOIN")


    if (earlystop == 1) {
      dselect[trial] <- mtd_est <-  NA
    }
    else {
      dselect[trial] <- mtd_est <- BOIN::select.mtd(target, n, y, cutoff.eli,
                                              extrasafe, offset)$MTD
    }
  }

  return(list(res_pat=res_pat,mtd_est=mtd_est))

}


# Helper function for 3+3

tpt_stopcheck <- function(res_pat,mtd_est,tox_rates) {

  prev_pat <- res_pat %>%
    dplyr::filter(dose_level==mtd_est)

  if (nrow(prev_pat) < 3) {

    # No patients have been treated at mtd_est -> continue at that dose level
    stop <- FALSE

  } else if (nrow(prev_pat) < 6) {

    if (sum(prev_pat$dlt) > 1) {

      # Dose de-escalation if more than 1 DLT
      mtd_est <- mtd_est - 1
      stop <- FALSE

    } else if (sum(prev_pat$dlt)==1) {

      # 3 Patients with 1 DLT -> include 3 more patients at that dose level
      stop <- FALSE


    } else {

      # 3 Patients with 0 DLT -> Dose escalation
      mtd_est <- mtd_est + 1
      stop <- FALSE
    }

  } else if (nrow(prev_pat) == 6) {

    if (sum(prev_pat$dlt) > 1) {

      # Dose de-escalation if more than 1 DLT
      mtd_est <- mtd_est - 1
      stop <- FALSE

    } else if (sum(prev_pat$dlt)<=1) {

      # 6 Patients with 1 or fewer DLT -> Dose escalation
      stop <- FALSE
      mtd_est <- mtd_est + 1

    }
  }




  # In case of previous dose de-escalation: stop trial if previous dose has been tested before
   prev_pat <- res_pat %>%
    dplyr::filter(dose_level==mtd_est)
  prev_pat_higher <- res_pat %>%
    dplyr::filter(dose_level>mtd_est)

  if (sum(prev_pat$dlt)>1) {
    stop <- TRUE
    mtd_est <- mtd_est - 1
  } else if (nrow(prev_pat)>0 & nrow(prev_pat_higher)>0) {
    stop <- TRUE
  }

  # Lowest dose is too toxic -> stop trial
  if (mtd_est==0) stop <- TRUE

  # Highest dose has been tested and is not too toxic <- stop trial
  if (mtd_est > length(tox_rates)) {
    stop <- TRUE
    mtd_est <- length(tox_rates)
  }

  return(list(mtd_est=mtd_est,stop=stop))
}


# Helper function for BCRM adapted from bcrm::bcrm

bcrm_ndose <- function(stop=list(nmax=NULL, nmtd=NULL, precision=NULL,
                                 nmin=NULL, safety=NULL),
                       data=NULL, p.tox0=NULL, sdose=NULL, dose=NULL, ff,
                       prior.alpha, cohort=3, target.tox, constrain=TRUE, only.below=FALSE,
                       sdose.calculate="mean", pointest="plugin",
                       tox.cutpoints=NULL, loss=NULL,
                       start=NULL, simulate=FALSE, nsims=1, truep=NULL,
                       threep3=FALSE, threep3.start = 1, threep3.esc.only=FALSE, method="exact", burnin.itr=2000,
                       production.itr=2000,
                       bugs.directory="c:/Program Files/WinBUGS14/",
                       plot=FALSE, seed=NULL,  quietly=10,  file=NULL,
                       N, tox, notox, quantiles = c(0.025, 0.25, 0.50, 0.75, 0.975)){


  # For compatability with older versions, allow logical quietly
  if (is.logical(quietly)){
    if (quietly) quietly <- nsims + 1
    else quietly <- 1
  }

  # Checks of argument inputs
  if(missing(N) & is.null(stop$nmax) & is.null(stop$nmtd) & is.null(stop$precision))
    stop("At least one stopping rule must be provided using the stop argument")
  if(!missing(N)){
    stop$nmax <- N
    warning("N is deprecated and users should now use the stop argument to specify the maximum sample size")
  }
  if(!missing(tox) | !missing(notox)){
    stop("tox and nontox arguments are deprecated and users should now use the data argument to specify previous data,  see ?bcrm")
  }
  if(!(length(stop$precision) %in% c(0, 2))) stop("stop$precision must be a vector of length two")
  if(!is.null(stop$nmax) & !is.null(stop$nmin)) {if(stop$nmin>stop$nmax) stop("stop$nmin must be less than stop$nmax")}
  if(!is.null(stop$safety)){
    if(stop$safety<=0 | stop$safety>=1) stop("stop$safety must be a probability between 0 and 1")
  }
  if(missing(p.tox0) & missing(sdose)) stop("Either p.tox0 or sdose must be specified")
  if(!missing(p.tox0) & !missing(sdose)) stop("Only one of p.tox0 and sdose must be specified")
  if(sdose.calculate!="mean" & sdose.calculate!="median") stop("sdose.calculate must be either `mean' or `median'")
  if((is.character(pointest) & pointest!="mean" & pointest!="plugin") | is.numeric(pointest) & (pointest<0 | pointest>1)) stop("pointest must be either `plugin',  `mean' or an EWOC feasibility quantile between 0 and 1")
  if(is.numeric(pointest) & method=="exact") stop("EWOC design must be fitted using MCMC methods")
  if(!is.null(tox.cutpoints) & method=="exact") stop("Escalation based on toxicity intervals must be fit using MCMC. Please specify either method=`rjags',  method='BRugs' or method='R2WinBUGS'")
  k<-max(length(p.tox0), length(sdose))
  if(simulate & is.null(truep)) stop("truep must be specified if simulating data")
  if(simulate){
    plot <- FALSE
  }
  if(!(method %in% c("exact", "rjags", "BRugs", "R2WinBUGS"))) stop("method must be either `exact',  `rjags',  `BRugs' or `R2WinBUGS'")
  ## Check to see if ff is one of "ht", "logit1", "power", "logit2"
  if((!ff %in% c("ht", "logit1", "power", "logit2"))) stop("ff must be one of `ht',  `logit1',  `power' or `logit2'")

  if(ff=="logit2" & method=="exact") warning("Exact method slow for 2-parameter model,  suggest using rjags (MCMC)")
  if(constrain & is.null(start) & is.null(data)) stop("A starting dose level must be specified using `start' if constrain==TRUE")
  if(!is.null(start)){if(constrain & is.null(data) & (start %in% 1:k)==FALSE) {stop("Starting dose required but not specified as one of 1, 2, ..., k.")}}
  if((!is.null(tox.cutpoints) & is.null(loss)) | (is.null(tox.cutpoints) & !is.null(loss))) stop("Both tox.cutpoints and loss must be specified to conduct escalation based on toxicity intervals")
  if(!is.null(tox.cutpoints) & length(loss)!=length(tox.cutpoints)+1) stop("The number of losses must be one more than the number of cutpoints")
  if(!is.null(tox.cutpoints)) pointest <- NULL
  if(!is.null(tox.cutpoints) & ff!="logit2") warning("One-parameter models are designed as working models only,  and should not be used with an escalation strategy based on intervals of the posterior probabilities of toxicity")

  if(ff=="logit2" & (length(prior.alpha[[2]])<2 | length(prior.alpha[[3]])<2)) stop("second and third components of `prior.alpha' must be vectors of size 2")

  if(threep3==TRUE & threep3.esc.only==TRUE & threep3.start!=1) stop("For 3+3 design with escalation only, starting dose for 3+3 must be equal to lowest dose (1)")
  if(threep3==TRUE){if(((!missing(sdose) & !(threep3.start %in% 1:length(sdose))) | (!missing(p.tox0) & !(threep3.start %in% 1:length(p.tox0))) )){stop("Start dose for 3+3 design must be one of the available dose levels")}}
  # Set seed if specified
  if(!is.null(seed)){
    set.seed(seed)
  }

  if(missing(sdose)){
    alpha.prior.plug <- if(prior.alpha[[1]]==1){
      ifelse(sdose.calculate=="mean", prior.alpha[[2]]*prior.alpha[[3]], median(getprior(prior.alpha,  10000)))
    } else if(prior.alpha[[1]]==2){
      0.5*(prior.alpha[[2]]+prior.alpha[[3]])
    } else if(prior.alpha[[1]]==3){
      ifelse(sdose.calculate=="mean", exp(prior.alpha[[2]]+prior.alpha[[3]]/2), exp(prior.alpha[[2]]))
    } else if(prior.alpha[[1]]==4){
      if(sdose.calculate=="mean"){exp(prior.alpha[[2]]+diag(prior.alpha[[3]])/2)} else {exp(prior.alpha[[2]])}
    }
    sdose <- bcrm::find.x(ff, p.tox0, alpha=alpha.prior.plug)
  }
  # Checking that length of truep in simulation study is same length as sdose
  if(simulate & length(truep)!=length(sdose)) stop("Length of truep must be the same as length of sdose or p.tox0")
  # Checking that length of dose (if specified) is same length as sdose
  if(length(dose)>0 & length(dose)!=length(sdose)) stop("Length of dose must be the same as length of sdose or p.tox0")
  # Check data contains the correct variables
  if(!is.null(data)){
    if(any(!(c("patient", "dose", "tox") %in% names(data)))) stop("data must have variables named 'patient',  'dose' and 'tox'")
    data <- data[order(data$patient), ]
    if(any(data$patient != 1:dim(data)[1])) stop("'patient' variable in data must be an ascending vector of positive integers")
    if(any(!(data$tox %in% c(0, 1)))) stop("'tox' variable in data must be a vector of zeros (no toxicity) and ones (toxicity)")
    if(any(!(data$dose %in% 1:length(sdose)))) stop(paste("'dose' variable in data must contain the dose levels (1 to ", length(sdose), ")", sep=""))
    #if(!is.null(start)) warning("start no longer needs to be specified if data is given; using last recruited patient as current dose")
    start <- as.numeric(data$dose[1])
  }

  ## Cannot calculate quantiles if method=="exact" so stop$precision and stop$safety cannot be used
  if(method=="exact" & (!is.null(stop$precision) | !is.null(stop$safety))){ stop("exact method cannot be used with stop$precision or stop$safety,  please use MCMC instead")}

  ## Allowing access to fast exact computation if method="exact" & simulate=TRUE & stopping rules do not depend on posterior quantiles
  if(method=="exact" & simulate & is.null(stop$precision)){ method <- "exact.sim" }
  ## If method=="exact" and a two-parameter model is fitted,  only relevant escalation posterior quantities are calculated (apart from trial end)
  if(method=="exact" & ff=="logit2"){
    method <- "exact.sim"
    if(plot){
      plot <- FALSE
      warning("Plot function not available for exact computation of 2-parameter model")
    }
  }

  k <- length(sdose)

  # Set up

  if(is.null(stop$safety)){
    quantiles <- quantiles
  } else {
    quantiles <- sort(unique(c(0.025, 0.25, 0.50, 0.75, 0.975, 1-stop$safety)))
  }


  if(is.null(data)){
    new.tox<- rep(0,k)
    new.notox<-rep(0,k)
    current<-start-1
    first<-TRUE
  }else{
    new.tox <- as.numeric(xtabs(tox ~ factor(dose, levels=1:k), data=data))
    new.notox <- as.numeric(xtabs((1 - tox) ~ factor(dose,levels=1:k), data=data))
    current <- as.numeric(data$dose[dim(data)[1]])
    first<-FALSE
  }

  alpha <-switch(method
                 ,rjags=bcrm::Posterior.rjags(new.tox,new.notox,sdose,ff,prior.alpha,burnin.itr,production.itr)
                 ,BRugs=bcrm::Posterior.BRugs(new.tox,new.notox,sdose,ff,prior.alpha,burnin.itr,production.itr)
                 ,R2WinBUGS=bcrm::Posterior.R2WinBUGS(new.tox,new.notox,sdose,ff,prior.alpha,burnin.itr,production.itr,bugs.directory)
                 ,exact=bcrm::Posterior.exact(new.tox,new.notox,sdose,ff,prior.alpha)
                 ,exact.sim=bcrm::Posterior.exact.sim(new.tox,new.notox,sdose,ff,prior.alpha,pointest))

  ncurrent <- sum(new.tox + new.notox)

  newdata <- data

  prior.ndose<-switch(method
                      ,rjags=bcrm:::nextdose(alpha,sdose,ff,target.tox,constrain,first,pointest,current,tox.cutpoints,loss,quantiles,only.below)
                      ,BRugs=bcrm:::nextdose(alpha,sdose,ff,target.tox,constrain,first,pointest,current,tox.cutpoints,loss,quantiles,only.below)
                      ,R2WinBUGS=bcrm:::nextdose(alpha,sdose,ff,target.tox,constrain,first,pointest,current,tox.cutpoints,loss,quantiles,only.below)
                      ,exact=bcrm:::nextdose.exact(alpha,sdose,ff,target.tox,constrain,first,pointest,current,only.below)
                      ,exact.sim=bcrm:::nextdose.exact.sim(alpha,sdose,ff,target.tox,constrain,first,pointest,current,only.below)
  )
  ndose <- prior.ndose

  stopped <- bcrm:::stop.check(stop, target.tox, ncurrent, ndose, new.tox, new.notox, simulate)
  ndose <- stopped$ndose

  results <- list(dose=dose, sdose=sdose, tox=new.tox, notox=new.notox, ndose=list(ndose), constrain=constrain, start=start, target.tox=target.tox, ff=ff, method=method, pointest=pointest, tox.cutpoints=tox.cutpoints, loss=loss, prior.alpha=prior.alpha, data=data)
  class(results) <- "bcrm"

  return(list(ndose=ndose$ndose,stop=stopped$stop))
}





