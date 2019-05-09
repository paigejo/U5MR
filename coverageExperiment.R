# this script is for calculating coverage for gaussian processes
library(fields)

set.seed(123)
n = 50 # number of observations
m = 100 # number of equally spaced prediction points at which we compute the coverage
sigma = 1 # GP marginal sd
phi = 1 / 10 # correlation scale
tau = 1 / 10 # nugget sd
kappa = 1 # matern smoothness
ss = seq(0, 3, l=m + 1)[-1] - 3 / (2*m) # s*, the prediction locations (33 is a last within the data domain)
ns = 100 # number of times the sample s
nu = 100 # number of times to sample u | s, s*
ny = 200 # number of times the sample y | u(s)

# precompute covariance matrix of the prediction locations
SigmaP = stationary.cov(matrix(ss, ncol=1), Covariance="Matern", theta=phi, smoothness=kappa)

# Function for compute conditional mean and variance for normal distribution given data.
# Xp | Xd has (conditional) mean:
# muc = muP + SigmaPtoD %*% SigmaD^(-1) %*% (Xd - muD)
# and (conditional) variance:
# Sigmac = SigmaP - SigmaPtoD %*% SigmaD^(-1) %*% SigmaDtoP
conditionalNormal0 = function(Xd, SigmaD, SigmaPtoD, SigmaDInv=NULL) {
  
  # SigmaDInv = solve(SigmaD) # NOTE: luckily we only need to do this once.
  # muc = muP + SigmaPtoD %*% SigmaDTildeInv %*% (Xd - muD)
  # Sigmac = SigmaP - SigmaPtoD %*% SigmaDInv %*% SigmaDtoP
  
  # compute conditional mean and variance of zeta
  # muc = muP + SigmaPtoD %*% solve(SigmaD, Xd - muD)
  if(is.null(SigmaDInv)) {
    muc = SigmaPtoD %*% solve(SigmaD, Xd) # mean is 0
    Sigmac = SigmaP - SigmaPtoD %*% solve(SigmaD, t(SigmaPtoD))
  }
  else {
    muc = SigmaPtoD %*% SigmaDInv %*% Xd # mean is 0
    Sigmac = SigmaP - SigmaPtoD %*% SigmaDInv %*% t(SigmaPtoD)
  }
  
  return(list(muc=muc, Sigmac=Sigmac))
}

sSamples = matrix(runif(n*ns), nrow=ns)
uSamples = array(dim=c(ns, n + m, nu))
uHats = array(dim=c(ns, m, nu, ny))
ySamples = array(dim=c(ns, n, nu, ny))
cov50 = array(dim=c(ns, nu, ny, m)) # all values are either 1 or 0
cov60 = array(dim=c(ns, nu, ny, m)) # all values are either 1 or 0
cov70 = array(dim=c(ns, nu, ny, m)) # all values are either 1 or 0
cov80 = array(dim=c(ns, nu, ny, m)) # all values are either 1 or 0
cov90 = array(dim=c(ns, nu, ny, m)) # all values are either 1 or 0
cov95 = array(dim=c(ns, nu, ny, m)) # all values are either 1 or 0
nnd = matrix(nrow=ns, ncol=m)
for(i in 1:ns) {
  print(paste0("number of s samples: ", i, "/", ns))
  
  # draw sample of s
  s = sSamples[i,]
  
  # calculate nearest neighbor distances for each prediction location
  nnd[i,] = apply(rdist(matrix(ss, ncol=1), matrix(s, ncol=1)), 1, function(x) {min(x[x != 0])})
  
  # generate covariance matrices, and invert data covariance
  SigmaPtoD = stationary.cov(matrix(ss, ncol=1), matrix(s, ncol=1), Covariance="Matern", theta=phi, smoothness=kappa)
  SigmaD = stationary.cov(matrix(s, ncol=1), Covariance="Matern", theta=phi, smoothness=kappa)
  diag(SigmaD) = diag(SigmaD) + rep(tau^2, n)
  Sigma = stationary.cov(matrix(c(s, ss), ncol=1), Covariance="Matern", theta=phi, smoothness=kappa)
  L = t(chol(Sigma))
  SigmaDInv = solve(SigmaD) # we're only inverting the matrix directly because it's small and we need to run this many times
  
  # generate samples of u | s, s*
  uSamples[i,,] = L %*% matrix(rnorm(nu * (n + m)), ncol=nu)
  
  for(j in 1:nu) {
    
    # draw sample of u | s, s*
    uSample = c(uSamples[i,,j])
    uss = uSample[-(1:n)]
    
    for(k in 1:ny) {
      
      # generate sample of y | u(s)
      ySample = uSample[1:n] + rnorm(n, sd=tau)
      ySamples[i,, j, k] = ySample
      
      # compute conditional distribution of u(s*) | y
      out = conditionalNormal0(ySample, SigmaD, SigmaPtoD, SigmaDInv)
      uHat = out$muc # predictions
      uHats[i,, j, k] = uHat
      sds = sqrt(diag(out$Sigmac)) # predictive standard deviations
      
      # determine if the true values of u(s*) are contained in the intervals
      cov50[i, j, k, ] = (uss >= uHat + qnorm(0.25, sd = sds)) & (uss <= uHat + qnorm(0.75, sd = sds))
      cov60[i, j, k, ] = (uss >= uHat + qnorm(0.2, sd = sds)) & (uss <= uHat + qnorm(0.8, sd = sds))
      cov70[i, j, k, ] = (uss >= uHat + qnorm(0.15, sd = sds)) & (uss <= uHat + qnorm(0.85, sd = sds))
      cov80[i, j, k, ] = (uss >= uHat + qnorm(0.1, sd = sds)) & (uss <= uHat + qnorm(0.9, sd = sds))
      cov90[i, j, k, ] = (uss >= uHat + qnorm(0.05, sd = sds)) & (uss <= uHat + qnorm(0.95, sd = sds))
      cov95[i, j, k, ] = (uss >= uHat + qnorm(0.025, sd = sds)) & (uss <= uHat + qnorm(0.975, sd = sds))
    }
  }
}

## aggregate in different ways
# average over y
cov50y = apply(cov50, c(1, 2, 4), mean)
cov60y = apply(cov60, c(1, 2, 4), mean)
cov70y = apply(cov70, c(1, 2, 4), mean)
cov80y = apply(cov80, c(1, 2, 4), mean)
cov90y = apply(cov90, c(1, 2, 4), mean)
cov95y = apply(cov95, c(1, 2, 4), mean)

# average over 

# save results
save(n, m, sigma, phi, tau, kappa, ss, ns, nu, ny, sSamples, uSamples, uHats, ySamples, 
     cov50, cov60, cov70, cov80, cov90, cov95, 
     file="coverageSimulation.RData")




