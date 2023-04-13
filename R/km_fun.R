km_fun <- function(t, time, ind){
  fit <- survfit(Surv(time, ind)~1)
  index <- sum(t >= fit$time) 
  prob <- fit$surv[index]
  return(prob)
}