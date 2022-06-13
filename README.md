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
seed <- 1
```

* **method** Vector of dose-finding designs to be employed. The designs will be applied sequentially with the estimated MTD from the previous design as the starting dose. Possible options are "Intra" for Intra-patient escalation, "TPT" for the 3+3 design, "BCRM" for (Bayesian) Continual Reassessment Method and "BOIN" for Bayesian Optimal Interval Design.
* **tox_rates** Vector of true toxicity rates for each dose level in ascending order
* **target** Target toxicity rate
* **intra_days** How many single treatments are given
* **nsim** Number of simulated trials when simulating multiple trials

### Simulating a single trial
``` r
set.seed(seed)
single_trial <- sim_trial(method=method,tox_rates=tox_rates,target=target,intra_days=intra_days)
single_trial$res_pat
# Result of simulated trial
# pat dose_level dlt design
# 1          1   0  Intra
# 1          2   0  Intra
# 1          3   0  Intra
# 1          4   1  Intra
# 1          3   0    TPT
# 2          3   1    TPT
# 3          3   1    TPT
# 4          2   0    TPT
# 5          2   0    TPT
# 6          2   0    TPT
# 7          2   0    TPT
# 8          2   0    TPT
# 9          2   0    TPT
```
* **pat** Patient number
* **dose_level** Dose level of the treatment
* **dlt** Whether a DLT occured (1=DLT, 0=no DLT)
* **design** Which design is used for the treatment

The first three rows display the intra-patient dose escalation stage for the first patient until they experience a DLT at dose level 4. Rows 4 to 13 display the results of the following stage which employs a 3+3 design. Patients 1, 2 and 3 start at dose level 3, i.e. the highest dose level with no DLT in the intra-patient dose escalation stage. As there are two DLTs at dose level 3 the dose is de-escalated to dose level 2. Patients 4, 5 and 6 don't experience any DLT and therefore dose level 2 is chosen as the MTD. A dose expansion phase follows and treats patients 7, 8 and 9 at dose level 2.
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
# [1] 2
#
# $res_dose
# # A tibble: 4 x 3
# dose_level     n   dlt
# <dbl> <dbl> <int>
# 1          1  0.25     0
# 2          2  6.25     0
# 3          3  3.25     2
# 4          4  0.25     1
```
* **n_pat** Displays the total number of patients included in the trial. In this case one patient in the intra-patient stage and nine patients in the 3+3 stage sums to ten. 
* **n_dlt** Displays the total number of DLTs in the trial. In this case one DLT in the intra-patient stage and two DLTs in the 3+3 stage sum to three.
* **mtd_est** Displays the estimated MTD at the end of the trial.
* **res_dose** Data frame which displays the number of fully treated patients (n) and the number of DLTS (dlt) at each dose level. A dose given during the intra-stage counts as 1/intra_days full treatments. In this case one patient in the intra-patient stage and six patients in the 3+3 stage sum to 6.25 full treatments at dose level 2.

### Simulating multiple trials
``` r
multiple_trials <- sim_trials(method=method,tox_rates=tox_rates,target=target,intra_days=intra_days,nsim=nsim,seed=seed)
multiple_trials$res
# Result of simulated trials
# n_pat n_dlt mtd_est
# 10     3       2
# 20     5       4
# 17     6       3
# 23     7       2
# 17     7       3
# 16     4       3
# 23     9       2
# 20     7       2
# 13     6       2
# 26    11       1
```
Each row displays the number of patients (n_pat), the number of DLTs (n_dlt) and the estimated MTD (mtd_est) for a simulated trial. Row 1 are the results of the example above of a single simulated trial.
``` r
summary_trials(multiple_trials)
# Summarized results of simulated trials
# $res_median
# n_pat n_dlt mtd_est
# apply.res_list.res..2..median.  18.5   6.5       2
```
res_median displays the median values for the number of patients (n_pat), the number of DLTs (n_dlt) and the estimated MTD (mtd_est) over the simulated trials.
``` r
# $accuracy
# [1] 0.3
```
accuracy is the percentage of trials which recommend the correct dose for phase II. 
``` r
# $freq_dose
# # A tibble: 5 x 4
# dose_level n_pat freq_pat true_mtd
# <dbl> <dbl>    <dbl> <lgl>   
#   1          1   8.5   0.0466 FALSE   
# 2          2  47.5   0.260  FALSE   
# 3          3  49.2   0.270  TRUE    
# 4          4  48.8   0.267  FALSE   
# 5          5  28.5   0.156  FALSE   
```
freq_dose displays the doses given at each dose level
* **dose_lvel** The dose level
* **n_pat** The total number of fully treated patients at each dose level. As before a dose in the intra-patient escalation stage is counted as 1/intra_days full treatments
* **freq_rp2d** The frequency a dose level has been given over all trials combined
* **true_mtd** Whether the dose level is the correct recommendation for phase II

``` r
# $freq_mtd_est
# # A tibble: 4 x 4
# mtd_est n_rp2d freq_rp2d true_mtd
# <dbl>  <int>     <dbl> <lgl>   
#   1       1      1       0.1 FALSE   
# 2       2      5       0.5 FALSE   
# 3       3      3       0.3 TRUE    
# 4       4      1       0.1 FALSE   
```
freq_mtd_est displays how often a dose is recommended for phase II
* **mtd_est** The dose level
* **n_rp2d** Total number of trials for which that dose level is recommendend for phase II
* **freq_rp2d** Percentage of trials for which that dose level is recommendend for phase II
* **true_mtd** Whether the dose level is the correct recommendation for phase II

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

