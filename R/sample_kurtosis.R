sample_kurtosis <- function(x){
  N = length(x)
  kur = N*(N+1)/((N-1)*(N-2)*(N-3))*sum((x-mean(x))^4)/var(x)^2 - 3*(N-1)^2/((N-2)*(N-3))
  
  return(kur)
}