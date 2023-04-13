survival.QC1 <- function(surv_time, ind, study_time){
  seed = 321 
  n = length(surv_time)
  
  # consider the exponential distribution
  rate_hat = sum(ind)/sum(surv_time)
  t.event = rexp(n, rate = rate_hat)
  t.ind = rep(1,times = n)
  exp_dis = weighted_distance(surv_time, ind, t.event, t.ind, study_time, 0, rate_hat)
  
  # consider the weibull distribution
  fit_weibull <- survreg(Surv(surv_time, ind)~1, dist = "weibull")
  shape_hat <- 1/fit_weibull$scale
  scale_hat <- exp(as.numeric(fit_weibull$coef[1]))
  set.seed(seed)
  t.event = rweibull(n, shape = shape_hat, scale = scale_hat)
  t.ind = rep(1,times = n)
  weibull_dis = weighted_distance(surv_time, ind, t.event, t.ind, study_time, 1, c(shape_hat, scale_hat))
  
  if (exp_dis < weibull_dis){
    cat(paste("The distribution of survival time is closer to Exponential distribution with rate =",round(rate_hat, digits = 3), 
              "compared to Weibull distribution with shape =", round(shape_hat, digits = 3), "and scale =", round(scale_hat, digits = 3), "\n"))
    
    out <- list(as.numeric(exp_dis), as.numeric(weibull_dis))
    names(out) <- c("exptest", "weibulltest")
    return(out)
    
  }
  else if (exp_dis==weibull_dis){
    cat(paste("The distribution of survival time has no difference between Exponential distribution with rate =",round(rate_hat, digits = 3), 
              "and Weibull distribution with shape =", round(shape_hat, digits = 3), "and scale =", round(scale_hat, digits = 3), "\n"))
    
    out <- list(as.numeric(exp_dis), as.numeric(weibull_dis))
    names(out) <- c("exptest", "weibulltest")
    return(out)
    
  }
  else{
    cat(paste("The distribution of survival time is closer to Weibull distribution with shape =", round(shape_hat, digits = 3), "and scale =", round(scale_hat, digits = 3), 
              "compared to Exponential distribution with rate =",round(rate_hat, digits = 3), "\n"))
    
    out <- list(as.numeric(exp_dis), as.numeric(weibull_dis))
    names(out) <- c("exptest", "weibulltest")
    return(out)
  }
}