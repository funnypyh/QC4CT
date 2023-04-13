# QC4CT

Title: Quality control for clinical trials 

Author: Zhanpeng Xu and Yeonhee Park

Maintainer: Yeonhee Park <ypark56@wisc.edu>

Description: This is an R package to check the trial assumptions such as accrual distribution and survival parametric distribution. There are three main functions. Function \code{accrual.QC} provides test statistics and p-value of Kolmogorov-Smirnov test for the identification of accrual distribution. Functions \code{survival.QC1} and \code{survival.QC2} provide weighted Kaplan-Meier statistics of Pepe and Fleming (1989) performing distance test for the identification of survival distribution. \code{survival.QC1} works when the observed survival time and censoring indicator are given; and \code{survival.QC2} works when the observed survival time and accrual time are given. 

Suggests: ExtDist, cubature, survival

License: GPL-2

Reference: Pepe, M. S. and Fleming, T. R. (1989). Weighted Kaplan-Meier statistics: a class of distance tests for censored survival data. Biometrics, 497-507.

## Example

library(ExtDist) 

y <- runif(100, min=0, max=1); z <- cumsum(y)

res1 <- accrual.QC(z)

y <- rgamma(100, shape=2, scale=0.5); z <- cumsum(y)

res2 <- accrual.QC(z)

library(cubature) 

library(survival) 

data(data1)

study_time <- max(data1[,1])+1

res <- survival.QC1(surv_time=data1[,1], ind=data1[,2], study_time)


data(data2)

res <- survival.QC2(surv_time=data2[,1], arr_time=data2[,2])

