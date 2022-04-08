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

multiple_trials <- sim_trials(method=method,tox_rates=tox_rates,target=target,intra_days=intra_day,nsim=nsim)
# Result of simulated trials
multiple_trials$res
# Summarized results of simulated trials
summary_trials(multiple_trials)
```

## Accessing simulated data

``` r
data(S1_BCRM)
summary_trials(S1_BCRM)
```

