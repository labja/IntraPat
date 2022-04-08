# IntraPat

## Installation

You can install this package from GitHub with
``` r
# install.packages("devtools")
devtools::install_github("labja/IntraPat")
```

## Example

``` r
library(IntraPat)

method <- c("Intra","TPT")
tox_rates <- seq(0.1,0.5,0.1)
target <- 0.3
intra_days <- 4
nsim <- 10

single_trial <- sim_trial(method=method,tox_rates=tox_rates,target=target,intra_days=intra_days)
# Result of simulated trial
single_trial$res_pat
# Summarized results of simulated trial
summary_trial(single_trial)

multiple_trials <- sim_trials(method=method,tox_rates=tox_rates,target=target,intra_days=intra_days,nsim=nsim)
# Result of simulated trials
multiple_trials$res
# Summarized results of simulated trials
summary_trials(multiple_trials)
```

## Accessing simulated data

``` r
data(S1_BCRM)
data(S1_BOIN)
data(S1_TPT)
data(S1_Intra_BCRM)
data(S1_Intra_BOIN)
data(S1_Intra_TPT)
data(S4_BCRM)
data(S4_BOIN)
data(S4_TPT)
data(S4_Intra_BCRM)
data(S4_Intra_BOIN)
data(S4_Intra_TPT)
data(S6_BCRM)
data(S6_BOIN)
data(S6_TPT)
data(S6_Intra_BCRM)
data(S6_Intra_BOIN)
data(S6_Intra_TPT)
data(S8_BCRM)
data(S8_BOIN)
data(S8_TPT)
data(S8_Intra_BCRM)
data(S8_Intra_BOIN)
data(S8_Intra_TPT)
```

