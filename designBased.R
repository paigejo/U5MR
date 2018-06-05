##################################################################################
# designBased.R                                                                  #
#    Produce county level estimates using pixel-wise spatial model with          #
#    aggregation using true population and BYM spatial model                     #
##################################################################################

# Load library
library(INLA)

# Load data
# load("simDataMulti.RData")
# load a different 1 of these depending on whether a cluster effect should be included 
# in the simulation of the data or not (tausq is the cluster effect variance)
load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOver2.RData")
# load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOver2.RData")

# include urban/rural effect
includeUrbanRural = TRUE

# Get true ratios of urban/rural
urbRatio = vector('numeric', length = 47)
for(i in 1:47){
    idx = which(as.numeric(factor(SRSDat$eaDat$admin1)) == i)
    urbRatio[i] = sum(SRSDat$eaDat$urban[idx])/length(idx)
}

# Define formula
if(includeUrbanRural) {
  formula = y ~ rural +
    f(idx, model="iid", 
      hyper=list(prec=list(param=c(0.5, 0.001488), prior="loggamma"))) + 
    f(idx2, model="besag",
      graph="Kenyaadm1.graph", 
      hyper=list(prec=list(param=c(0.5, 0.00360), prior="loggamma"))) +
    f(idxEps, model = "iid",
      hyper = list(prec = list(prior = 'loggamma', param = c(1,0.01))))
} else {
  formula = y ~ f(idx, model="iid", 
                  hyper=list(prec=list(param=c(0.5, 0.001488), prior="loggamma"))) + 
    f(idx2, model="besag",
      graph="Kenyaadm1.graph", 
      hyper=list(prec=list(param=c(0.5, 0.00360), prior="loggamma"))) +
    f(idxEps, model = "iid",
      hyper = list(prec = list(prior = 'loggamma', param = c(1,0.01))))
}


# Number of simulations for producing results
Nsim = 1000

# Help functions
logit = function(x){
    return(log(x/(1-x)))
}
expit = function(x){
    return(1/(1+exp(-x)))
}

# Go through datasets for SRSDat
sampCountySRSDat = array(NA, dim = c(47, 1000, 100))
for(i in 1:100){
    # Extract data
    currData = SRSDat$clustDat[[i]]
    currData$admin1 = factor(currData$admin1)
    
    # INLA data
    dat = list(y = currData$died,
               Ntrials = currData$numChildren,
               rural = 1-currData$urban,
               idx = as.numeric(currData$admin1),
               idx2 = as.numeric(currData$admin1),
               idxEps = 1:length(currData$died))
    
    # Add unobserved data to make sampling easier
    dat$y = c(rep(NA, 47*2), dat$y)
    dat$Ntrials = c(rep(25, 47*2), dat$Ntrials)
    dat$rural = c(rep(c(0,1), each = 47), dat$rural)
    dat$idx = c(rep(1:47, 2), dat$idx)
    dat$idx2 = c(rep(1:47, 2), dat$idx2)
    dat$idxEps = c(rep(NA, 47*2), dat$idxEps)
  
    # Run model
    result = inla(formula = formula,
                  family="binomial",
                  Ntrials = Ntrials,
                  data=dat,
                  control.compute = list(config = TRUE))
    
    # Simulate from posterior
    if(includeUrbanRural) {
      samp = inla.posterior.sample(n = Nsim, result = result)
      sampRural = matrix(NA, nrow = 47, ncol = Nsim)
      sampUrban = matrix(NA, nrow = 47, ncol = Nsim)
      for(j in 1:Nsim){
        sampRural[, j] = samp[[j]]$latent[47 + (1:47)]
        sampUrban[, j] = samp[[j]]$latent[1:47]
      }
      sampCountySRSDat[,,i] = logit(expit(sampUrban)*urbRatio + expit(sampRural)*(1-urbRatio))
    } else {
      samp = inla.posterior.sample(n = Nsim, result = result)
      result = matrix(NA, nrow = 47, ncol = Nsim)
      for(j in 1:Nsim){
        result[, j] = samp[[j]]$latent[1:47]
      }
      sampCountySRSDat[,,i] = result
    }
}

# Go through datasets for overSampDat
sampCountyOverSampDat = array(NA, dim = c(47, 1000, 100))
for(i in 1:100){
    # Extract data
    currData = overSampDat$clustDat[[i]]
    currData$admin1 = factor(currData$admin1)
    
    # INLA data
    dat = list(y = currData$died,
               Ntrials = currData$numChildren,
               rural = 1-currData$urban,
               idx = as.numeric(currData$admin1),
               idx2 = as.numeric(currData$admin1),
               idxEps = 1:length(currData$died))
    
    # Add unobserved data to make sampling easier
    dat$y = c(rep(NA, 47*2), dat$y)
    dat$Ntrials = c(rep(25, 47*2), dat$Ntrials)
    dat$rural = c(rep(c(0,1), each = 47), dat$rural)
    dat$idx = c(rep(1:47, 2), dat$idx)
    dat$idx2 = c(rep(1:47, 2), dat$idx2)
    dat$idxEps = c(rep(NA, 47*2), dat$idxEps)
    
    # Run model
    result = inla(formula = formula,
                  family="binomial",
                  Ntrials = Ntrials,
                  data=dat,
                  control.compute = list(config = TRUE))
    
    # Simulate from posterior
    if(includeUrbanRural) {
      samp = inla.posterior.sample(n = Nsim, result = result)
      sampRural = matrix(NA, nrow = 47, ncol = Nsim)
      sampUrban = matrix(NA, nrow = 47, ncol = Nsim)
      for(j in 1:Nsim){
        sampRural[, j] = samp[[j]]$latent[47 + (1:47)]
        sampUrban[, j] = samp[[j]]$latent[1:47]
      }
      sampCountyOverSampDat[,,i] = logit(expit(sampUrban)*urbRatio + expit(sampRural)*(1-urbRatio))
    } else {
      samp = inla.posterior.sample(n = Nsim, result = result)
      result = matrix(NA, nrow = 47, ncol = Nsim)
      for(j in 1:Nsim){
        result[, j] = samp[[j]]$latent[1:47]
      }
      sampCountyOverSampDat[,,i] = result
    }
}

processSamples = function(samp){
    # 80% credible intervals
    CI = t(apply(X = samp, MARGIN = 1, FUN = quantile, probs = c(0.1, 0.5, 0.9)))
    mm = rowMeans(samp)
    ss = apply(X = samp, MARGIN = 1, FUN = sd)
    
    return(list(logit = list(CI = CI,
                             mean = mm,
                             stddev = ss),
                prob = list(CI = exp(CI))))
}

## SRSdata
    Q10 = matrix(NA, nrow = 47, ncol = 100)
    Q50 = matrix(NA, nrow = 47, ncol = 100)
    Q90 = matrix(NA, nrow = 47, ncol = 100)
    mm = matrix(NA, nrow = 47, ncol = 100)
    ss = matrix(NA, nrow = 47, ncol = 100)
    for(i in 1:100){
        tmp = processSamples(sampCountySRSDat[,,i])
        Q10[,i] = tmp$logit$CI[,1]
        Q50[,i] = tmp$logit$CI[,2]
        Q90[,i] = tmp$logit$CI[,3]
        mm[,i] = tmp$logit$mean
        ss[,i] = tmp$logit$stddev
    }
    resSRSdat = list(Q10 = Q10,
                     Q50 = Q50,
                     Q90 = Q90,
                     mean = mm,
                     stddev = ss)
    
## overSampDat
    Q10 = matrix(NA, nrow = 47, ncol = 100)
    Q50 = matrix(NA, nrow = 47, ncol = 100)
    Q90 = matrix(NA, nrow = 47, ncol = 100)
    mm = matrix(NA, nrow = 47, ncol = 100)
    ss = matrix(NA, nrow = 47, ncol = 100)
    for(i in 1:100){
        tmp = processSamples(sampCountyOverSampDat[,,i])
        Q10[,i] = tmp$logit$CI[,1]
        Q50[,i] = tmp$logit$CI[,2]
        Q90[,i] = tmp$logit$CI[,3]
        mm[,i] = tmp$logit$mean
        ss[,i] = tmp$logit$stddev
    }
    resOverSRSdat = list(Q10 = Q10,
                         Q50 = Q50,
                         Q90 = Q90,
                         mean = mm,
                         stddev = ss)

# Full result
    designRes = list(SRSdat = resSRSdat,
                     overSampDat = resOverSRSdat)
    # save(file = 'kenyaSpatialDesignResultNew.RData', designRes = designRes)
    save(file = paste0('kenyaSpatialDesignResultNewTausq0UrbRur', 
                         includeUrbanRural, '.RData'), designRes = designRes)
    # save(file = paste0('kenyaSpatialDesignResultNewTausq0.01UrbRur', 
    #                    includeUrbanRural, '.RData'), designRes = designRes)
    