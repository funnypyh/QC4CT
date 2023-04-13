accrual.QC <- function(x){
  
  time = diff(x)
  n = length(time)
  lnx = log(time)
  
  theta_hat = (sum(time*lnx)*n - sum(time)*sum(lnx))/(n*n)
  k_hat = mean(time)/theta_hat
  theta = n/(n-1)*theta_hat
  k = k_hat - (3*k_hat-2/3*k_hat/(1+k_hat)-4/5*k_hat/(1+k_hat)^2)/n
  
  gam_test = ks.test(time, "pgamma", shape = k, scale = theta)
  gam_dis = gam_test$statistic
  gam_pval = gam_test$p.value
  
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
  
  para1 = min(time) 
  para2 = max(time)
  uni_test = ks.test(time, "punif", para1, para2)
  uni_dis = uni_test$statistic
  uni_pval = uni_test$p.value
  
  beta_test = ks.test(time, "pBeta_ab", shape1 = shape1_hat, shape2 = shape2_hat,
                      a = a_hat, b= b_hat)
  beta_dis = beta_test$statistic
  beta_pval = beta_test$p.value
  
  para = mean(time)
  poi_test = ks.test(time, "pexp", 1/para)
  poi_dis = poi_test$statistic
  poi_pval = poi_test$p.value
  
  gammatest <- c(as.numeric(gam_dis), gam_pval)
  exptest <- c(as.numeric(poi_dis), poi_pval)
  betatest <- c(as.numeric(beta_dis), beta_pval)
  uniformtest <- c(as.numeric(uni_dis), uni_pval)
  out <- list(gammatest, exptest, betatest, uniformtest)
  names(out) <- c("gammatest", "exptest", "betatest", "unifromtest")
  return(out)
  
}