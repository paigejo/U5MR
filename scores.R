# truth: realised value
# my.est: estimator
# my.var: variance
# lower: lower CI limit
# upper: upper CI limit
library(PearsonDS)
logit <- function(x) {
  log(x/(1-x))
}

expit <- function(x) {
  exp(x)/(1+exp(x))
}

mse <- function(truth, my.est, logit=TRUE, my.var=NULL, nsim=10, n=1){
  
  if(!logit){
    # this is the MSE of the median prediction
    # truth = expit(truth)
    # my.est = expit(my.est)
    
    # calculate the MSE of the mean prediction:
    # first draw different values of p from the predictive distribution
    if(is.null(my.var))
      stop("variance must be supplied when calculating mse on non-logit scale")
    quantiles = matrix(rep(seq(0, 1, l=nsim + 2)[c(-1, -(nsim + 2))], length(my.est)), ncol=nsim, byrow=TRUE)
    logitP = matrix(qnorm(quantiles, my.est, sd=sqrt(my.var)), ncol=nsim)
    my.est = expit(logitP)
    res = mean(sweep(my.est, 1, truth, "-")^2 * n^2)
  }
  else {
    res <- mean((truth - my.est)^2)
  }
  
  return(res)
}

bias <- function(truth, my.est, logit=TRUE, my.var=NULL, n=1){
  
  if(!logit) {
    # this is the bias of the median prediction, not the mean
    # truth = expit(truth)
    # my.est = expit(my.est)
    
    # compute the bias of the mean prediction
    if(is.null(my.var))
      stop("logit variance must be included for bias calculation on non-logit scale")
    my.est = logitNormMean(cbind(my.est, sqrt(my.var)))
  }
  res <- my.est - truth
  
  if(!logit) {
    return(mean(n * res))
  }
  
  return(mean(res))
}

crpsNormal <- function(truth, my.est, my.var, logit=TRUE, nsim=10, resultType="county", n=25){
  if(resultType == "EA" || resultType == "pixel")
    logit=FALSE
  
  if(logit) {
    sig = sqrt(my.var)
    x0 <- (truth - my.est) / sig
    res <- sig * (1 / sqrt(pi) -  2 * dnorm(x0) - x0 * (2 * pnorm(x0) - 1))
    
    ## sign as in Held (2008)
    res <- -res
  }
  else {
    # first draw different values of p from the predictive distribution
    quantiles = matrix(rep(seq(0, 1, l=nsim + 2)[c(-1, -(nsim + 2))], length(my.est)), ncol=nsim, byrow=TRUE)
    logitP = matrix(qnorm(quantiles, my.est, sd=sqrt(my.var)), ncol=nsim)
    p = expit(logitP)
    
    # for each value of p, calculate the crps and return the mean
    crpsVals = crpsBinomial(truth, p, n=n)
    res = mean(crpsVals)
  }
  
  res
}

# for a fixed binomial, determine the CRPS. 
# trueProportion: can be a vector
# p.est: can be a matrix, with number of rows equal to length(trueProportion)
# n: can be a vector, with length equal to length(trueProportion)
# empiricalProportion: if true, calculate the CRPS of the empirical proportion distribution on [0, 1] 
#                      rather than the binomial distribution on [0, n]
crpsBinomial <- function(trueProportion, p.est, n=25, parClust=NULL, empiricalProportion=TRUE) {
  trueCounts = trueProportion * n
  p.est = matrix(p.est, nrow=length(trueProportion))
  
  # make sure n is a vector
  if(length(n) == 1)
    n = rep(n, length(trueProportion))
  
  rowFun = function(rowI) {
    toN = eval(0:n[rowI])
    thisCount = trueCounts[rowI]
    
    crpsBinomialSingle = function(colI) {
      thisP = p.est[rowI, colI]
      
      sum((pbinom(toN, n[rowI], thisP) - (thisCount <= toN))^2)
    }
    
    lapply(1:ncol(p.est), crpsBinomialSingle)
  }
  
  if(is.null(parClust))
    crpsByRegion = lapply(1:nrow(p.est), function(i) {mean(unlist(rowFun(i)))})
  else
    crpsByRegion = parLapply(parClust, 1:nrow(p.est), function(i) {mean(unlist(rowFun(i)))})
  
  if(empiricalProportion)
    sum(unlist(crpsByRegion) / n)
  else
    sum(unlist(crpsByRegion))
}

# for a discrete random variable, determine the CRPS. 
# probs: a vector of probabilities of length n + 1
crpsCounts <- function(trueCount, probs, parClust=NULL, empiricalProportion=TRUE) {
  n = length(probs) - 1
  toN = eval(0:n)
  
  # get cdf
  cdf = cumsum(probs)
  
  if(is.null(parClust))
    res = sum((cdf - (trueCount <= toN))^2)
  else {
    res = parSapply(parClust, toN, function(x) {trueCount <= x})
    res = parSapply(parClust, 1:(n + 1), function(x) {(cdf[x] - res[x])^2})
    res = sum(res)
  }
  
  if(empiricalProportion)
    res / n
  else
    res
}

dss = function(truth, my.est, my.var){
  
  my.sd = sqrt(my.var)
  # proportional to logarithmic score assuming normality
  dss_score <- (2*log(my.sd) + ((truth - my.est)/my.sd)^2)
  
  return(mean(dss_score))
}

coverage = function(truth, lower, upper, logit = TRUE){
  if(!logit){
    lower = expit(lower)
    upper = expit(upper)
    # truth = expit(truth)
  }
  
  if(any(lower >= upper)) {
    warning("lower >= upper, reordering")
    tmp = lower
    wrongOrder = lower >= upper
    lower[wrongOrder] = upper[wrongOrder]
    upper[wrongOrder] = tmp[wrongOrder]
  }
  
  res = mean(lower <= truth & upper >= truth)
  
  return(res)
}

intervalWidth = function(lower, upper, logit = TRUE) {
  if(!logit){
    lower = expit(lower)
    upper = expit(upper)
  }
  
  res = mean(abs(upper - lower))
  if(any(lower > upper))
    warning("some interval widths are negative, taking absolute value")
  res
}

# taken from logitnorm package.  Calculates the mean of a distribution whose 
# logit is Gaussian.
logitNormMean = function(muSigmaMat, parClust=NULL, ...) {
  if(length(muSigmaMat) > 2) {
    if(is.null(parClust))
      apply(muSigmaMat, 1, logitNormMean, ...)
    else
      parApply(parClust, muSigmaMat, 1, logitNormMean, ...)
  }
  else {
    mu = muSigmaMat[1]
    sigma = muSigmaMat[2]
    fExp <- function(x) plogis(x) * dnorm(x, mean = mu, sd = sigma)
    integrate(fExp, -Inf, Inf, abs.tol = 0, ...)$value
  }
}

# calculates the E[UV] for U and V bernoulli random variables with the given distribution
# of p on the logit scale
logitNormCrossExpectation = function(muSigmaMat=matrix(c(0, 1), ncol=2), ...) {
  
  # now calculate the expectation of the product of two of the random variables
  if(length(muSigmaMat) > 2) {
    if(is.null(parClust))
      apply(muSigmaMat, 1, logitNormMean, ...)
    else
      parApply(parClust, muSigmaMat, 1, logitNormMean, ...)
  }
  else {
    mu = muSigmaMat[1]
    sigma = muSigmaMat[2]
    f11 <- function(x) plogis(x)^2 * dnorm(x, mean = mu, sd = sigma)
    p11 = integrate(f11, -Inf, Inf, abs.tol = 0, ...)$value
    p11
  }
}

# calculates the cor(U, V) for U and V bernoulli random variables with the given distribution
# of p on the logit scale
logitNormCor = function(muSigmaMat=matrix(c(0, 1), ncol=2)) {
  # first calculate marginal expectations
  expectedProbs = logitNormMean(muSigmaMat)
  
  # calculate the cross expectations
  crossExpectations = logitNormCrossExpectation(muSigmaMat)
  
  # now calculate the correlations
  numerator = crossExpectations - expectedProbs^2
  denominator = expectedProbs * (1 - expectedProbs)
  numerator / denominator
}

# generate a binomial confidence interval with at least the significance requested. The
# interval is generated by starting at the maximum probability count and expanding 
# outward to the highest probability counts.
# NOTE: the return confidence interval is on counts, not the index of the probs vector
generateBinomialInterval = function(probs, significance = .80) {
  if(any(is.na(probs)))
    return(list(lower=NA, upper=NA, leftRejectProb=NA, rightRejectProb=NA))
  # Find the maximum probability, end expand outward from there
  maxI = which.max(probs)
  leftI = maxI
  rightI = maxI
  
  expandInterval = function(leftI, rightI) {
    if(length(leftI) > 1 || length(rightI) > 1) {
      print("woah there")
    }
    currentProb = sum(probs[leftI:rightI])
    
    if(currentProb >= significance)
      return(list(leftI = leftI, rightI = rightI))
    
    if(leftI == 1)
      leftProb = 0
    else
      leftProb = probs[leftI - 1]
    if(rightI == length(probs))
      rightProb = 0
    else
      rightProb = probs[rightI + 1]
    
    if(rightProb == leftProb) {
      # in this case, expand symmetrically to the mode, or randomly if distances from mode are the same
      modeI = which.max(probs)
      leftDist = modeI - leftI
      rightDist = rightI - modeI
      if(leftDist == rightDist) {
        if(runif(1) >= .5)
          return(list(leftI = leftI, rightI = rightI + 1))
        else
          return(list(leftI = leftI - 1, rightI = rightI))
      } else if(leftDist > rightDist)
        return(list(leftI = leftI, rightI = rightI + 1))
      else
        return(list(leftI = leftI - 1, rightI = rightI))
    } else if(rightProb >= leftProb)
      return(list(leftI = leftI, rightI = rightI + 1))
    else
      return(list(leftI = leftI - 1, rightI = rightI))
  }
  
  # keep expanding the interval until were done (i doesn't matter here)
  for(i in 1:length(probs)) {
    out = expandInterval(leftI, rightI)
    if(leftI == out$leftI && rightI == out$rightI)
      break
    
    leftI = out$leftI
    rightI = out$rightI
  }
  
  # determine probability of rejection at boundary
  totalProb = sum(probs[leftI:rightI])
  extraProb = totalProb - significance
  leftProb = probs[leftI]
  rightProb = probs[rightI]
  if(leftProb <= rightProb) {
    leftRejectProb = extraProb
    rightRejectProb = 0
  }
  else {
    rightRejectProb = extraProb
    leftRejectProb = 0
  }
  
  # return results (subtract one to the index to get the count)
  c(lower=leftI - 1, upper=rightI - 1, leftRejectProb=leftRejectProb, rightRejectProb=rightRejectProb)
}

# approximate the sum of binomial distributions using the method from:
# https://stackoverflow.com/questions/15926448/approximate-the-distribution-of-a-sum-of-binomial-random-variables-in-r
# http://www.dtic.mil/dtic/tr/fulltext/u2/a266969.pdf
# use the first 4 moments of the sum, fits with a pearson distribution
dSumBinom = function(k, ns=25, ps=.5) {
  if(length(ns) == 1)
    ns = rep(ns, length(ps))
  if(length(ps) == 1)
    ps = rep(ps, length(ns))
  
  # get the first 4 cumulants
  k.1<-sum(ns*ps)
  k.2<-sum(ns*ps*(1-ps))
  k.3<-sum(ns*ps*(1-ps)*(1-2*ps))
  k.4<-sum(ns*ps*(1-ps)*(1-6*ps*(1-ps)))
  
  # obtain skewness and excess kurtosis
  beta.1<-k.3^2/k.2^3
  beta.2<-k.4/k.2^2
  
  # get the moments, and return the probability mass
  # (by integrating the density of the continuous pearson distribution)
  moments <- c(mean=k.1,variance=k.2,skewness=sqrt(beta.1),kurtosis=beta.2 + 3)
  if(!identical(k, 0:sum(ns))) {
    results = k
    zeros = k == 0
    results[zeros] = ppearson(.5,moments=moments)
    results[!zeros] = ppearson(k[!zeros] + .5,moments=moments) - ppearson(k[!zeros] - .5,moments=moments)
    results
  } else {
    # the commented out code is technically correct, but the code below is a good approximation that saves time
    # results = ppearson(k + .5,moments=moments)
    # results[2:length(results)] = diff(results)
    # results
    results = dpearson(k, moments=moments)
    results * (1 / sum(results))
  }
}

# approximate the sum of binomials where the probabilities or sampled from a random distribution
# using a Pearson distribution
# pMat: a matrix of probabilities, the number of rows equaling the length of ns, and the number 
#       of columns corresponding to the number of samples of p for each discrete distribution
dSumBinomRandom = function(k, ns=25, pMat, parClust=NULL) {
  if(is.null(parClust) || length(ns) <= 100 || ncol(pMat) <= 30)
    dSumBinoms = matrix(apply(pMat, 2, dSumBinom, k=k, ns=ns), ncol=ncol(pMat))
  else
    dSumBinoms = matrix(parApply(parClust, pMat, 2, dSumBinom, k=k, ns=ns), ncol=ncol(pMat))
  rowMeans(dSumBinoms)
}

# approximate the sum of binomials where the probabilities or sampled from a random distribution 
# using simulations from binomial distributions
# pMat: a matrix of probabilities, the number of rows equaling the length of ns, and the number 
#       of columns corresponding to the number of samples of p for each discrete distribution
dSumBinomRandomSim = function(k, ns=25, pMat, nSim = 1000) {
  dSumBinoms = matrix(apply(pMat, 2, dSumBinomSim, k=k, ns=ns, nSim=nSim), ncol=ncol(pMat))
  rowMeans(dSumBinoms)
}

# same as dSumBinom except we used simulations from binomial distributions conditional on p 
# to approximate the distribution rather than using the pearson distribution approximation
dSumBinomSim = function(k, ns=25, ps=.5, nSim=1000) {
  if(length(ns) == 1)
    ns = rep(ns, length(ps))
  if(length(ps) == 1)
    ps = rep(ps, length(ns))
  
  countSims = matrix(rbinom(length(ns) * nSim, ns, ps), nrow=length(ns), ncol=nSim)
  cdf = ecdf(colSums(countSims))
  cdf(k) - cdf(k - 1)
}


