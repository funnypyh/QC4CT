sample_skewness <- function(x){
  N = length(x)
  skew = N/(N-1)/(N-2)*sum((x-mean(x))^3)/sqrt(var(x)^3)
  
  return(skew)
}