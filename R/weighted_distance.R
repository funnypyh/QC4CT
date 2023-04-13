weighted_distance <- function(surv1, ind1, surv2, ind2, study_time, acc_dis, acc_para){
  
  n1 = length(surv1)
  n2 = length(surv2)
  n = n1+n2
  
  p1 = n1/n
  p2 = n2/n
  c_ind1 = 1-ind1
  c_ind2 = 1-ind2
  
  weight_fun <- function(t){
    mole = km_fun(t, surv1, c_ind1)*km_fun(t, surv2, c_ind2)
    deno = p1*km_fun(t, surv1, c_ind1)+p2*km_fun(t, surv2, c_ind2)
    return(mole/deno)
  }
  
  int1 <- function(t){
    if (acc_dis==0){ # exponential distribution
      return(weight_fun(t)*(km_fun(t,surv1,ind1) - 1 + pexp(t, rate = acc_para[1])))
    }
    if (acc_dis==1){ # weibull distribution
      return(weight_fun(t)*(km_fun(t,surv1,ind1) - 1 + pweibull(t, shape = acc_para[1], scale = acc_para[2])))
    }
    
  }
  dis = sqrt(n1*n2/n)*adaptIntegrate(int1, 0, study_time)$integral 
  return(abs(dis))
}