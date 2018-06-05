# truth: realised value
# my.est: estimator
# my.var: variance
# lower: lower CI limit
# upper: upper CI limit

logit <- function(x){
  log(x/(1-x))
}

expit <- function(x){
  exp(x)/(1+exp(x))
}

mse <- function(truth, my.est, logit=TRUE){

  if(!logit){
    truth = expit(truth)
    my.est = expit(my.est)
  }
  
  res <- (truth - my.est)^2
  return(res)
}

bias <- function(truth, my.est, logit=TRUE){
  
  if(!logit){
    truth = expit(truth)
    my.est = expit(my.est)
  }
  res <- (truth-my.est)
  
  return(res)
}

crpsNormal <- function(truth, my.est, my.var){
  
  sig = sqrt(my.var)
  x0 <- (truth - my.est) / sig
  res <- sig * (1 / sqrt(pi) -  2 * dnorm(x0) - x0 * (2 * pnorm(x0) - 1))
  
  ## sign as in Held (2008)
  res <- -res
  
  return(res)
}

dss = function(truth, my.est, my.var){
  
  my.sd = sqrt(my.var)
  # proportional to logarithmic score assuming normality
  dss_score <- (2*log(my.sd) + ((truth - my.est)/my.sd)^2)
  
  return(dss_score)
}

coverage = function(truth, lower, upper, logit=TRUE){
  
  if(!logit){
    lower = expit(lower)
    upper = expit(upper)
    truth = expit(truth)
  }
  
  res = mean(lower <= truth & upper >= truth)
  
  return(res)
}
