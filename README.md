# IntraPat

## Installation

You can install this package from GitHub with
``` r
# install.packages("devtools")
devtools::install_github("labja/IntraPat")
```

## Example

### Setting parameters 
``` r
library(IntraPat)
method <- c("Intra","TPT")
tox_rates <- seq(0.1,0.5,0.1)
target <- 0.3
intra_days <- 4
nsim <- 10
seed <- 12345
```
Alternative design choices are "BCRM" and "BOIN".
### Simulating a single trial
``` r
single_trial <- sim_trial(method=method,tox_rates=tox_rates,target=target,intra_days=intra_days,seed=seed)
single_trial$res_pat
# Result of simulated trial
# pat dose_level dlt design
#   1          1   0  Intra
#   1          2   0  Intra
#   1          3   1  Intra
#   1          2   1    TPT
#   2          2   1    TPT
#   3          2   0    TPT
#   4          1   0    TPT
#   5          1   0    TPT
#   6          1   0    TPT
#   7          1   0    TPT
#   8          1   0    TPT
#   9          1   0    TPT
```
The first three rows display the intra-patient dose escalation stage for the first patient until they experience a DLT at dose level 3. Rows 4 to 12 display the results of the following stage which employs a 3+3 design. Patients 1, 2 and 3 start at dose level 2, i.e. the highest dose level with no DLT in the intra-patient dose escalation stage. As there are two DLTs at dose level 2 the dose is de-escalated to dose level 1. Patients 4, 5 and 6 don't experience any DLT and therefore dose level is chosen as the MTD. A dose expansion phase follows and treats patients 7, 8 and 9 at dose level 1.
``` r
# Summarized results of simulated trial
summary_trial(single_trial)
# $n_pat
# [1] 10
# 
# $n_dlt
# [1] 3
# 
# $mtd_est
# [1] 1
# 
# $res_dose
# # A tibble: 3 x 3
# dose_level     n   dlt
# <dbl> <dbl> <int>
# 1          1  6.25     0
# 2          2  3.25     2
# 3          3  0.25     1
```
n_pat displays the total number of patients included in the trial. In this case one patient in the intra-patient stage and nine patients in the 3+3 stage sums to ten. n_dlt displays the total number of DLTs in the trial. In this case one DLT in the intra-patient stage and two DLTs in the 3+3 stage sum to three. mtd_est displays the estimated MTD at the end of the trial. res_dose is a data frame which display the number of patients fully treated at each dose level (n) and the number of DLTS (dlt). A dose given during the intra-stage counts as 1/intra_days full treatments. In this case one patient in the intra-patient stage and six patients in the 3+3 stage sum to 6.25 full treatments at dose level 1.

``` r
multiple_trials <- sim_trials(method=method,tox_rates=tox_rates,target=target,intra_days=intra_days,nsim=nsim,seed=seed)
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

