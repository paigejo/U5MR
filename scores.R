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

# this function computes all scoring rules, with and without binomial variation
# truth: the true empirical proportions of mortality rates within the regions or enumeration areas of interest
# numChildren: a vector containing the number of children with the same length as truth
# logitEst: the predictive mean on the logit scale, of the same length as truth
# logitVar: the predictive variance on the logit scale, of the same length as truth
# est: the estimate on the probability scale, of the same length as truth. By default calculated from probMat
# var: the estimate variance on the probability scale, of the same length as truth. By default calculated from probMat
# logitL: the lower end of the credible interval on the logit scale
# logitU: the upper end of the credible interval on the logit scale
# probMat: a matrix of joint draws of probability values, with number of rows equal to the length of truth, a number of 
#          columns equal to the number of draws. If not included, a gaussian distribution is assumed on the logit scale.
# significance: the significance level of the credible interval. By default 80%
# nsim: the number of sample quantiles from the gaussian distribution when approximating the proportion scale 
#       distributions from the logit scale
# NOTE: Discrete, count level credible intervals are estimated based on the input probMat along with coverage and CRPS
getScores = function(truth, numChildren, logitEst, logitVar, est=NULL, var=NULL, probMat=NULL, significance=.8, nsim=10, 
                     do1row=TRUE) {
  # first calculate bias (the same with and without binomial variation). Since on empirical proportion scale, set n to 1
  thisBias = bias(truth, logitEst, logit=FALSE, logitVar, n=1)
  
  # calculate average predictive variance on the logit scale
  # thisLogitVar = mean(logitVar)
  
  ##### The below code was commented out because it is done automatically in the coverage function:
  # # if necessary, calculate credible interval boundaries on the logit scale without binomial variation
  # if(is.null(logitL))
  #   logitL = qnorm((1 - significance) / 2, logitEst, sqrt(logitVar))
  # if(is.null(logitU))
  #   logitU = qnorm(1 - (1 - significance) / 2, logitEst, sqrt(logitVar))
  # 
  # # if necessary, do the same with binomial variation based on the joint simulations in the probability matrix
  # if(is.null(probL))
  #   probL = apply(probMat, 1, function(x) {quantile(x, prob=(1 - significance) / 2)})
  # if(is.null(probU))
  #   probU = apply(probMat, 1, function(x) {quantile(x, prob=1 - (1 - significance) / 2)})
  # 
  # # calculate credible interval width on the proportion scale with and without binomial variation
  # thisWidthNoBinom = mean(expit(logitU) - expit(logitL))
  # thisWidthBinom = mean(probU - probL)
  
  # calculate mean squared error at the proportion scale
  if(is.null(est))
    thisMSE = mse(truth, logitEst, logit=FALSE, logitVar, nsim=nsim)
  else
    thisMSE = mse(truth, logit(est), logit=FALSE, logitVar, nsim=nsim)
  
  # calculate predictive variance at the proportion scale based on the joint simulation probability matrix 
  # if available or assume gaussian predictive distribution on the logit scale otherwise
  if(!is.null(probMat)) {
    if(is.null(var))
      var = apply(probMat, 1, sd)^2
  }
  else {
    quantiles = matrix(rep(seq(0, 1 - 1 / nsim, by = 1/nsim) + 1 / (2 * nsim), length(logitEst)), ncol=nsim, byrow=TRUE)
    logitP = matrix(qnorm(quantiles, logitEst, sd=sqrt(logitVar)), ncol=nsim)
    pMat = expit(logitP)
    if(is.null(var))
      var = apply(pMat, 1, sd)^2
  }
  
  # calculate predictive variance on the proportion scale
  thisVar = mean(var)
  
  # calculate coverage and credible interval width with and without binomial variation
  coverageNoBinom = coverage(truth, doLogit=FALSE, bVar=FALSE, numChildren=numChildren, logitEst=logitEst, logitVar=logitVar, 
                             probMat=probMat, significance=significance, returnIntervalWidth=TRUE, nsim=nsim)
  coverageBinom = coverage(truth, doLogit=FALSE, bVar=TRUE, numChildren=numChildren, logitEst=logitEst, logitVar=logitVar, 
                           probMat=probMat, significance=significance, returnIntervalWidth=TRUE, returnPMFs=TRUE, nsim=nsim)
  
  thisCoverageNoBinom = coverageNoBinom[1]
  thisCoverageBinom = coverageBinom$intervalInfo[1]
  thisWidthNoBinom = coverageNoBinom[2]
  thisWidthBinom = coverageBinom$intervalInfo[2]
  
  # also get the probability mass functions for the binomial predictive distributions for the CRPS calculations
  pmfs = coverageBinom$pmfs
  
  # calculate CRPS on the proportion scale with and without binomial variation
  thisCRPSNoBinom = crpsNormal(truth, logitEst, logitVar, logit=FALSE, n=numChildren, bVar=FALSE, probMat=probMat, nsim=nsim)
  thisCRPSBinom = crpsNormal(truth, logitEst, logitVar, logit=FALSE, n=numChildren, bVar=TRUE, probMat=probMat, nsim=nsim, pmfs=pmfs)
  
  # return the results
  if(do1row) {
    results = matrix(c(thisBias, thisVar, thisMSE, thisCRPSNoBinom, thisCRPSBinom, thisCoverageNoBinom, thisCoverageBinom, 
                       thisWidthNoBinom, thisWidthBinom), nrow=1)
    colnames(results) = c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")
    as.data.frame(results)
  } else {
    results = matrix(c(thisBias, thisVar, thisMSE, thisCRPSNoBinom, thisCoverageNoBinom, thisWidthNoBinom, 
                       thisBias, thisVar, thisMSE, thisCRPSBinom, thisCoverageBinom, thisWidthBinom), 
                     nrow=2, byrow = TRUE)
    colnames(results) = c("bias", "var", "mse", "crps", "coverage", "length")
    as.data.frame(results)
  }
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
    quantiles = matrix(rep(seq(0, 1 - 1 / nsim, by = 1/nsim) + 1 / (2 * nsim), length(my.est)), ncol=nsim, byrow=TRUE)
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

# truth: a vector of observations on the desired scale
# my.est: a vector of logit-scale predictions of the same length as truth 
# my.var: a vector of logit-scale predictive variances of the same length as truth
# logit: whether CRPS should be computed on the logit scale or the empirical proportion scale
# nsim: if binomial variation is included, the number of probabilities over which to integrate in the pmf approximation
# n: the binomial sample size. Can be a vector of the same length as truth
# bVar: whether or not binomial variation should be included in the predictive distribution
# probMat: if included, use these probability draws in the pmf integration/approximation. If the estimates are normal 
#          on the logit scale, do not include this since just using the estimates and variance will be accurate.
# pmfs: a list of probability mass functions of the same length as truth
crpsNormal <- function(truth, my.est, my.var, logit=TRUE, nsim=10, n=25, bVar=TRUE, probMat=NULL, pmfs=NULL){
  # if(resultType == "EA" || resultType == "pixel")
  #   logit=FALSE
  
  if(logit) {
    sig = sqrt(my.var)
    x0 <- (truth - my.est) / sig
    res <- sig * (1 / sqrt(pi) -  2 * dnorm(x0) - x0 * (2 * pnorm(x0) - 1))
    
    ## sign as in Held (2008)
    res <- -res
  }
  else {
    if(bVar) {
      # first draw different values of p from the predictive distribution if necessary
      if(is.null(probMat) && is.null(pmfs)) {
        quantiles = matrix(rep(seq(0, 1 - 1 / nsim, by = 1/nsim) + 1 / (2 * nsim), length(my.est)), ncol=nsim, byrow=TRUE)
        logitP = matrix(qnorm(quantiles, my.est, sd=sqrt(my.var)), ncol=nsim)
        probMat = expit(logitP)
      }
      
      # for each value of p, calculate the crps and return the mean
      crpsVals = crpsBinomial(truth, probMat, n=n, pmfs=pmfs)
      res = mean(crpsVals)
    } else {
      # no binomial variation is included. Integrate numerically on the proportion scale
      
      # compute the crps for this row of truth
      crpsRow = function(rowI) {
        thisTruth = truth[rowI]
        
        # either build the predictive cdf assuming normality on the logit scale or from the 
        # empirical distribution given by probMat if the user supplies it
        if(is.null(probMat)) {
          thisEst = my.est[rowI]
          thisVar = my.var[rowI]
          thisCdf = function(ws) {pnorm(logit(ws), thisEst, sqrt(thisVar))}
        } else {
          thisCdf = ecdf(probMat[rowI,])
        }
        
        intFun = function(ws) {
          (thisCdf(ws) - (ws >= thisTruth))^2
        }
        
        if(is.null(probMat)) {
          # when integrating we will set bounds on the integral to be reasonable to avoid 
          # faulty "divergent integral" error. The bounds will be 20 standard errors out 
          # of the estimate, making sure to include the truth.
          lowerBound = max(0, min(thisTruth - .01, expit(thisEst - 20 * sqrt(thisVar))))
          upperBound = min(1, max(thisTruth + .01, expit(thisEst + 20 * sqrt(thisVar))))
          integrate(intFun, lowerBound, upperBound)$value
        }
        else {
          # since we are using the empirical distribution, there is a closed form for the integral
          ps = probMat[rowI,]
          allPoints = sort(c(ps, thisTruth))
          deltas = diff(allPoints)
          sum(deltas * intFun(allPoints[1:length(ps)]))
        }
      }
      crpsVals = sapply(1:length(truth), crpsRow)
      
      res = mean(crpsVals)
    }
  }
  
  res
}

# for a set of independent fixed binomials, determine the CRPS.
# trueProportion: can be a vector
# p.est: can be a matrix, with number of rows equal to length(trueProportion). In that case, the 
#        number of columns can be thought of as the number of draws of p in the case of a random p.
# n: can be a vector, with length equal to length(trueProportion)
# empiricalProportion: if true, calculate the CRPS of the empirical proportion distribution on [0, 1] 
#                      rather than the binomial distribution on [0, n]
# pmfs: a list of probability mass functions with the same length as trueProportion
crpsBinomial <- function(trueProportion, p.est=NULL, n=25, parClust=NULL, empiricalProportion=TRUE, pmfs=NULL) {
  # get observations at the count scale, make sure n is a vector
  trueCounts = trueProportion * n
  if(length(n) == 1)
    n = rep(n, length(trueProportion))
  
  if(is.null(pmfs)) {
    p.est = matrix(p.est, nrow=length(trueProportion))
  }
  else
    cdfs = lapply(pmfs, function(pmf) {cumsum(pmf)})
  
  rowFun = function(rowI) {
    thisCount = eval(trueCounts[rowI])
    toN = eval(0:n[rowI])
    
    if(is.null(pmfs)) {
      ## in this case we must compute the probability mass function for this region ourselves
      thisN = eval(n[rowI])
      
      crpsBinomialSingle = function(colI) {
        thisP = p.est[rowI, colI]
        sum((pbinom(toN, thisN, thisP) - (thisCount <= toN))^2)
      }
      
      # return the region CRPS, taking the expectation over the distribution of p
      mean(unlist(sapply(1:ncol(p.est), crpsBinomialSingle)))
    } else
      sum((cdfs[[rowI]] - (thisCount <= toN))^2)
  }
  
  if(is.null(parClust))
    crpsByRegion = lapply(1:length(trueProportion), function(i) {rowFun(i)})
  else
    crpsByRegion = parLapply(parClust, 1:length(trueProportion), function(i) {rowFun(i)})
  
  if(empiricalProportion)
    mean(unlist(crpsByRegion) / n) # multiply by delta p for the discrete integral
  else
    mean(unlist(crpsByRegion))
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

# for a discrete random variable, determine the CRPS. 
# probs: a vector of probabilities of length n + 1
crpsCountsNoBVar <- function(trueProportion, logitEst, logitVar, parClust=NULL, empiricalProportion=TRUE) {
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

# either include both lower and upper, or include either: 
#    - the probability draw matrix
#    - estimates and variances on the logit scale
# lower: the logit scale lower end of the credible interval
# upper: the logit scale upper end of the credible interval
# nsim: the number of draws from the distribution of p in the case that the logit estimates and variance are provided
#       to estimate the probability mass function
coverage = function(truth, lower=NULL, upper=NULL, doLogit = TRUE, bVar=FALSE, numChildren=NULL, probMat=NULL, logitEst=NULL, logitVar=NULL, 
                    nsim=10, significance=.8, returnIntervalWidth=FALSE, returnPMFs=FALSE){
  
  if(any(is.null(lower)) || any(is.null(upper))) {
    # if the user did not supply their own credible intervals, we must get them ourselves given the other information
    
    if(is.null(probMat) && (is.null(logitEst) || is.null(logitVar)))
      stop("either include both lower and upper, probMat, or both logitEst and logitVar")
    
    if(!is.null(logitEst) && !is.null(logitVar) && bVar) {
      # in this case, we must generate the probabilities over which we will later integrate to get the 
      # pmfs and credible intervals
      
      ## first we estimate the probability mass function
      # get the probabilities over which we will integrate to get the probability mass function
      # quantiles = matrix(rep(seq(0, 1, l=nsim + 2)[c(-1, -(nsim + 2))], length(logitEst)), ncol=nsim, byrow=TRUE)
      quantiles = matrix(rep(seq(0, 1 - 1 / nsim, by = 1/nsim) + 1 / (2 * nsim), length(logitEst)), ncol=nsim, byrow=TRUE)
      logitP = matrix(qnorm(quantiles, logitEst, sd=sqrt(logitVar)), ncol=nsim)
      probMat = expit(logitP)
    }
    else if(!is.null(logitEst) && !is.null(logitVar) && !bVar) {
      # We don't want to include binomial variation, and we have the gaussian distribution parameters on the logit scale. 
      # Use those parameters to generate the credible intervals
      
      lower = qnorm((1 - significance) / 2, logitEst, sqrt(logitVar))
      upper = qnorm(1 - (1 - significance) / 2, logitEst, sqrt(logitVar))
    }
    else {
      # we don't have information about the predictive distribution on the logit scale, and don't assume normality. 
      # Instead, use the user supplied to probability matrix probMat
      
      if(!bVar) {
        # we aren't accounting for binomial variation, so just take the quantiles of the probability draws
        CIs = logit(apply(probMat, 1, function(ps) {quantile(ps, probs=c((1 - significance) / 2, 1 - (1 - significance) / 2))}))
        lower = CIs[1,]
        upper = CIs[2,]
      }
      else {
        # don't need to do anything here, since it will be already done to probMat in the following if statement
      }
    }
    
    if(bVar) {
      # for each region, take the average binomial probability mass, averaging over that drawn probabilities (integrate)
      # do this for all possible count values to get the probability mass for each count
      getPMF = function(rowI) {
        ps = probMat[rowI,]
        n = numChildren[rowI]
        
        rowMeans(sapply(ps, function(p) {dbinom(0:n, n, p)}))
      }
      pmfs = lapply(1:nrow(probMat), getPMF)
      
      ## now that we have the pmfs, we can calculate the credible intervals (on the logit rather than count scale)
      # NOTE: it is fine if we have a 0 or 1 interval boundary, since if we take the logit and get infinite results, 
      #       comparisons still work
      # intervals = sapply(pmfs, getHDI, significance=significance)
      intervals = sapply(pmfs, getQuantileInterval, significance=significance)
      lower = logit(intervals[1,] / numChildren)
      upper = logit(intervals[2,] / numChildren)
      lowerRejectProb = intervals[3,]
      upperRejectProb = intervals[4,]
    }
  }
  
  # if(!doLogit){
  lower = expit(lower)
  upper = expit(upper)
  # truth = expit(truth)
  # }
  
  if(any(lower >= upper)) {
    warning("lower >= upper, reordering")
    tmp = lower
    wrongOrder = lower >= upper
    lower[wrongOrder] = upper[wrongOrder]
    upper[wrongOrder] = tmp[wrongOrder]
  }
  
  res = lower <= truth & upper >= truth
  
  if(returnIntervalWidth)
    width = upper - lower
  
  if(bVar) {
    # in the case of discrete credible intervals, reject if the truth is at the boundary with some probability
    atLower = truth == lower
    atUpper = truth == upper
    rejectLower = runif(sum(atLower)) < lowerRejectProb[atLower]
    rejectUpper = runif(sum(atUpper)) < upperRejectProb[atUpper]
    res[atLower] = res[atLower] - rejectLower
    res[atUpper] = res[atUpper] - rejectUpper
    
    # adjust interval width for random rejection if necessary
    if(returnIntervalWidth) {
      width[atLower] = width[atLower] - rejectLower / numChildren[atLower]
      width[atUpper] = width[atUpper] - rejectUpper / numChildren[atUpper]
    }
  }
  
  allResults = c(coverage=mean(res))
  if(returnIntervalWidth)
    allResults = c(allResults, width=mean(width))
  if(returnPMFs)
    allResults = list(intervalInfo=allResults, pmfs=pmfs)
  
  allResults
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
# logit is Gaussian. Each row of muSigmaMat has a mean and standard deviation 
# on the logit scale
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

# generate a highest density interval with at least the significance requested
# NOTE: the returned confidence interval is on counts, not the index of the probs vector
getHDI = function(probs, significance = .80) {
  if(any(is.na(probs)))
    return(list(lower=NA, upper=NA, leftRejectProb=NA, rightRejectProb=NA))
  
  # Find the maximum probability, and expand outward from there
  maxI = which.max(probs)
  n = length(probs) - 1
  cumulativeProbs = cumsum(probs)
  
  getIntervalWidth = function(leftI) {
    # totalProbs = probs[leftI] + cumulativeProbs[leftI:(n + 1)] - cumulativeProbs[leftI]
    # width = binarySearchMatch(1, totalProbs >= significance)
    width = binarySearchMatch(1, cumulativeProbs[leftI:(n + 1)], function(x, i) {(x[i] - x[1] + probs[leftI]) >= significance})
    rightI = leftI - 1 + width
    c(width = width, leftI = leftI, rightI = rightI, prob=cumulativeProbs[rightI] - cumulativeProbs[leftI] + probs[leftI])
  }
  
  # calculate interval widths for relevant starting points
  maxRelevantI = binarySearchMatch(1, cumulativeProbs >= 1 - significance)
  intervalMatrix = sapply(1:maxRelevantI, getIntervalWidth)
  
  # determine which interval is the highest density interval, and get its information
  # NOTE: there may be multiple highest density intervals, since the distribution is discrete. 
  #       Taking any of them is fine.
  minI = which.min(intervalMatrix[1,])
  leftI = intervalMatrix[2, minI]
  rightI = intervalMatrix[3, minI]
  totalProb = intervalMatrix[4, minI]
  
  # determine probability of rejection at boundary
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

# generate a credible interval with at least the significance requested by including utmost (1 - significance) / 2
# probability mass on each side
# NOTE: the returned confidence interval is on counts, not the index of the probs vector
getQuantileInterval = function(probs, significance = .80) {
  if(any(is.na(probs)))
    return(list(lower=NA, upper=NA, leftRejectProb=NA, rightRejectProb=NA))
  
  # get the cdf
  n = length(probs) - 1
  cumulativeProbs = cumsum(probs)
  
  # get left and right endpoints
  leftI = max(binarySearchMatch(1, cumulativeProbs >= (1 - significance) / 2) - 1, 1)
  probLeft = cumulativeProbs[leftI] - probs[leftI]
  extraProbLeft = (1 - significance) / 2 - probLeft
  rightI = max(binarySearchMatch(1, cumulativeProbs >= 1 - (1 - significance) / 2), 1)
  probRight = 1 - cumulativeProbs[leftI]
  extraProbRight = (1 - significance) / 2 - probRight
  
  # determine probability of rejection at boundaries
  leftRejectProb = extraProbLeft
  rightRejectProb = extraProbRight
  
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



