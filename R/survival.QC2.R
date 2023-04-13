survival.QC2 <- function(surv_time, arr_time){
  
  # surv_time = Survival time
  # arr_time = Accrual time
  
  acc_quality_control <- function(x){
    sample_skewness <- function(x){
      N = length(x)
      skew = N/(N-1)/(N-2)*sum((x-mean(x))^3)/sqrt(var(x)^3)
      
      return(skew)
    }
    sample_kurtosis <- function(x){
      N = length(x)
      kur = N*(N+1)/((N-1)*(N-2)*(N-3))*sum((x-mean(x))^4)/var(x)^2 - 3*(N-1)^2/((N-2)*(N-3))
      
      return(kur)
    }
    
    time = diff(x)
    n = length(time)
    lnx = log(time)
    
    theta_hat = (sum(time*lnx)*n - sum(time)*sum(lnx))/(n*n)
    k_hat = mean(time)/theta_hat
    theta = n/(n-1)*theta_hat
    k = k_hat - (3*k_hat-2/3*k_hat/(1+k_hat)-4/5*k_hat/(1+k_hat)^2)/n
    
    gam_test = ks.test(time, "pgamma", shape = k, scale = theta)
    gam_dis = gam_test$statistic
    
    # check the assumption for the beta distribution, mark = 0 means the assumption holds
    mark = 1
    skw = sample_skewness(time)
    kur = sample_kurtosis(time)
    if (skw == 0){
      if ( -2 < kur & kur < 0 ){
        mark = 0
        mu_hat = 3*(kur+2)/(-kur)
        shape1_hat = (3/2*kur+3)/(-kur)
        shape2_hat = (3/2*kur+3)/(-kur)
        
        ba_diff = sqrt(var(time))/2*sqrt((2+mu_hat)^2*skw^2+16*(1+mu_hat))
        a_hat = mean(time) - shape1_hat/mu_hat*ba_diff
        b_hat = a_hat+ba_diff
      }
    }
    else{
      if (skw^2-2 < kur & kur < 3/2*skw^2 ){
        mark = 0
        mu_hat = 3*(kur - skw^2+2)/(3/2*skw^2-kur)
        if (skw > 0){
          shape1_hat = mu_hat/2*(1-1/(sqrt(1+(16*(mu_hat+1))/(mu_hat+2)^2/skw^2)))
          shape2_hat = mu_hat/2*(1+1/(sqrt(1+(16*(mu_hat+1))/(mu_hat+2)^2/skw^2)))
        }
        else{
          shape1_hat = mu_hat/2*(1+1/(sqrt(1+(16*(mu_hat+1))/(mu_hat+2)^2/skw^2)))
          shape2_hat = mu_hat/2*(1-1/(sqrt(1+(16*(mu_hat+1))/(mu_hat+2)^2/skw^2)))
        }
        
        ba_diff = sqrt(var(time))/2*sqrt((2+mu_hat)^2*skw^2+16*(1+mu_hat))
        a_hat = mean(time) - shape1_hat/mu_hat*ba_diff
        b_hat = a_hat+ba_diff
      }
    }
    
    if (mark == 1){
      para1 = min(time) 
      para2 = max(time)
      uni_test = ks.test(time, "punif", para1, para2)
      uni_dis = uni_test$statistic
    }
    else{
      beta_test = ks.test(time, "pBeta_ab", shape1 = shape1_hat, shape2 = shape2_hat,
                          a = a_hat, b= b_hat)
      beta_dis = beta_test$statistic
    }
    
    if ( mark == 1){
      if (gam_dis < uni_dis){
        return(c(0, gam_dis, uni_dis, mark, k, theta, para1, para2))
      }
      else if (gam_dis==uni_dis){
        return(c(0, gam_dis, uni_dis, mark, k, theta, para1, para2)) # we only use the gamma process
      }
      else{
        return(c(1, gam_dis, uni_dis, mark, k, theta, para1, para2))
      }
    }
    
    if ( mark == 0){
      if (gam_dis < beta_dis){
        return(c(0, gam_dis, beta_dis, mark, k, theta, shape1_hat, shape2_hat, a_hat, b_hat))
      }
      else if (gam_dis==beta_dis){
        return(c(0, gam_dis, beta_dis, mark, k, theta, shape1_hat, shape2_hat, a_hat, b_hat)) # we only use the gamma process
      }
      else{
        return(c(2, gam_dis, beta_dis, mark, k, theta, shape1_hat, shape2_hat, a_hat, b_hat))
      }
    }
  }
  
  seed = 321
  n = length(surv_time)
  arrival.t <- arr_time 
  t.event <- surv_time
  tobs = arrival.t[n]+1
  ind=t.ind <- ifelse(arrival.t+t.event == tobs,0,1)
  accdis <- acc_quality_control(arrival.t)
  acc_dis <- accdis[1]
  acc_para <- accdis[-c(1, 2, 3, 4)]
  
  # accrual time
  set.seed(seed)
  if (acc_dis==0){ # gamma distribution
    wait.t = rgamma(n,shape = acc_para[1], scale = acc_para[2])
    arrival.t = cumsum(wait.t)
  }
  if (acc_dis==1){ # uniform distribution
    wait.t = runif(n,min = acc_para[3], max = acc_para[4])
    arrival.t = cumsum(wait.t)
  }
  if (acc_dis==2){ # beta distribution
    wait.t = rBeta_ab(n,shape1 = acc_para[3], shape2 = acc_para[4], a = acc_para[5], b = acc_para[6])
    arrival.t = cumsum(wait.t)
  }
  
  # consider the exponential distribution
  rate_hat = sum(ind)/sum(surv_time)
  event.t = rexp(n, rate = rate_hat)
  tobs = arrival.t[n]
  t.event = ifelse(arrival.t[1:n]+
                     event.t[1:n]<=tobs,
                   event.t[1:n]
                   ,tobs-arrival.t[1:n])
  t.ind = ifelse(arrival.t[1:n]+
                   event.t[1:n]<=tobs,1,0)
  study_time = tobs 
  exp_dis = weighted_distance(surv_time, ind, t.event, t.ind, study_time, 0, rate_hat)
  
  # consider the weibull distribution
  fit_weibull <- survreg(Surv(surv_time, ind)~1, dist = "weibull")
  shape_hat <- 1/fit_weibull$scale
  scale_hat <- exp(as.numeric(fit_weibull$coef[1]))
  event.t = rweibull(n, shape = shape_hat, scale = scale_hat)
  tobs = arrival.t[n]
  t.event = ifelse(arrival.t[1:n]+
                     event.t[1:n]<=tobs,
                   event.t[1:n]
                   ,tobs-arrival.t[1:n])
  t.ind = ifelse(arrival.t[1:n]+
                   event.t[1:n]<=tobs,1,0)
  study_time = max(c(tobs, arr_time[n]))
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