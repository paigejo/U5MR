##################################################################################
# designBased.R                                                                  #
#    Produce county level estimates using pixel-wise spatial model with          #
#    aggregation using true population and BYM spatial model                     #
##################################################################################

runBYM = function(tausq=0.1^2, test=FALSE, includeUrbanRural=TRUE, includeCluster=TRUE) {
  
  # load and relevant data
  if(!test)
    load(paste0("simDataMultiBeta-1.75margVar0.0225tausq", round(tausq, 4), "gamma-1HHoldVar0urbanOverSamplefrac0.RData"))
  else
    load(paste0("simDataMultiBeta-1.75margVar0.0225tausq", round(tausq, 4), "gamma-1HHoldVar0urbanOverSamplefrac0Test.RData"))
  # load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOverSamplefrac0.25.RData")
  # load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOverSamplefrac0.25Test.RData")
  # load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOverSamplefrac0.25Test.RData")
  
  # Get true ratios of urban/rural
  urbRatio = vector('numeric', length = 47)
  for(i in 1:47){
    idx = which(as.numeric(factor(SRSDat$eaDat$admin1)) == i)
    urbanI = idx[SRSDat$eaDat$urban[idx]]
    urbRatio[i] = sum(SRSDat$eaDat$numChildren[urbanI])/sum(SRSDat$eaDat$numChildren[idx])
  }
  
  # Define formula
  if(includeUrbanRural) {
    if(includeCluster) {
      formula = y ~ rural +
        f(idx, model="iid", 
          hyper=list(prec=list(param=c(0.5, 0.001488), prior="loggamma"))) + 
        f(idx2, model="besag",
          graph="Kenyaadm1.graph", 
          hyper=list(prec=list(param=c(0.5, 0.00360), prior="loggamma"))) +
        f(idxEps, model = "iid",
          hyper = list(prec = list(prior = 'loggamma', param = c(1,0.01))))
    } else {
      formula = y ~ rural + 
        f(idx, model="iid", 
          hyper=list(prec=list(param=c(0.5, 0.001488), prior="loggamma"))) + 
        f(idx2, model="besag",
          graph="Kenyaadm1.graph", 
          hyper=list(prec=list(param=c(0.5, 0.00360), prior="loggamma")))
    }
  } else {
    if(includeCluster) {
      formula = y ~ f(idx, model="iid", 
                      hyper=list(prec=list(param=c(0.5, 0.001488), prior="loggamma"))) + 
        f(idx2, model="besag", graph="Kenyaadm1.graph", 
          hyper=list(prec=list(param=c(0.5, 0.00360), prior="loggamma"))) +
        f(idxEps, model = "iid",
          hyper = list(prec = list(prior = 'loggamma', param = c(1,0.01))))
    } else {
      formula = y ~ 
        f(idx, model="iid", 
          hyper=list(prec=list(param=c(0.5, 0.001488), prior="loggamma"))) + 
        f(idx2, model="besag", graph="Kenyaadm1.graph", 
          hyper=list(prec=list(param=c(0.5, 0.00360), prior="loggamma")))
    }
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
  sampCountySRSDat = array(NA, dim = c(47, Nsim, length(SRSDat$clustDat)))
  sampCountySRSDatMod = array(NA, dim = c(47, Nsim, length(SRSDat$clustDat)))
  for(i in 1:length(SRSDat$clustDat)){
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
    dat$Ntrials = c(rep(1, 47*2), dat$Ntrials)
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
    
    if(includeCluster)
      sampClusterSigmaSRS = 1:Nsim
    
    # Simulate from posterior
    if(includeUrbanRural) {
      samp = inla.posterior.sample(n = Nsim, result = result)
      sampRural = matrix(NA, nrow = 47, ncol = Nsim)
      sampUrban = matrix(NA, nrow = 47, ncol = Nsim)
      sampRuralMod = matrix(NA, nrow = 47, ncol = Nsim)
      sampUrbanMod = matrix(NA, nrow = 47, ncol = Nsim)
      
      for(j in 1:Nsim){
        sampRural[, j] = samp[[j]]$latent[47 + (1:47)]
        sampUrban[, j] = samp[[j]]$latent[1:47]
        if(includeCluster) {
          # if cluster effect is included, must debias predictions in each modeled strata
          clusterSigma = sqrt(1/samp[[j]]$hyperpar[3])
          muSigmaMatRural = cbind(sampRural[, j], clusterSigma)
          muSigmaMatUrban = cbind(sampUrban[, j], clusterSigma)
          sampRuralMod[, j] = logitNormMean(muSigmaMat = muSigmaMatRural)
          sampUrbanMod[, j] = logitNormMean(muSigmaMat = muSigmaMatUrban)
        }
      }
      sampCountySRSDat[,,i] = logit(expit(sampUrban)*urbRatio + expit(sampRural)*(1-urbRatio))
      sampCountySRSDatMod[,,i] = logit(sampUrbanMod*urbRatio + sampRuralMod*(1-urbRatio))
    } else {
      samp = inla.posterior.sample(n = Nsim, result = result)
      sampCounty = matrix(NA, nrow = 47, ncol = Nsim)
      sampCountyMod = matrix(NA, nrow = 47, ncol = Nsim)
      for(j in 1:Nsim){
        sampCounty[, j] = samp[[j]]$latent[1:47]
        if(includeCluster) {
          # if cluster effect is included, must debias predictions in each modeled strata
          clusterSigma = sqrt(1/samp[[j]]$hyperpar[3])
          muSigmaMat = cbind(sampCounty[, j], clusterSigma)
          sampCountyMod[, j] = logitNormMean(muSigmaMat = muSigmaMat)
        }
      }
      sampCountySRSDat[,,i] = sampCounty
      sampCountySRSDatMod[,,i] = logit(sampCountyMod)
    }
  }
  
  # Go through datasets for overSampDat
  sampCountyOverSampDat = array(NA, dim = c(47, Nsim, length(overSampDat$clustDat)))
  sampCountyOverSampDatMod = array(NA, dim = c(47, Nsim, length(overSampDat$clustDat)))
  for(i in 1:length(overSampDat$clustDat)) {
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
    dat$Ntrials = c(rep(1, 47*2), dat$Ntrials)
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
    
    if(includeCluster) {
      sampClusterSigmaOverSampDat = 1:Nsim
    }
    
    # Simulate from posterior
    if(includeUrbanRural) {
      samp = inla.posterior.sample(n = Nsim, result = result)
      sampRural = matrix(NA, nrow = 47, ncol = Nsim)
      sampUrban = matrix(NA, nrow = 47, ncol = Nsim)
      sampRuralMod = matrix(NA, nrow = 47, ncol = Nsim)
      sampUrbanMod = matrix(NA, nrow = 47, ncol = Nsim)
      for(j in 1:Nsim){
        sampRural[, j] = samp[[j]]$latent[47 + (1:47)]
        sampUrban[, j] = samp[[j]]$latent[1:47]
        
        if(includeCluster) {
          # if cluster effect is included, must debias predictions in each modeled strata
          clusterSigma = sqrt(1/samp[[j]]$hyperpar[3])
          muSigmaMatRural = cbind(sampRural[, j], clusterSigma)
          muSigmaMatUrban = cbind(sampUrban[, j], clusterSigma)
          sampRuralMod[, j] = logitNormMean(muSigmaMat = muSigmaMatRural)
          sampUrbanMod[, j] = logitNormMean(muSigmaMat = muSigmaMatUrban)
        }
      }
      sampCountyOverSampDat[,,i] = logit(expit(sampUrban)*urbRatio + expit(sampRural)*(1-urbRatio))
      sampCountyOverSampDatMod[,,i] = logit(sampUrbanMod*urbRatio + sampRuralMod*(1-urbRatio))
    } else {
      samp = inla.posterior.sample(n = Nsim, result = result)
      sampCounty = matrix(NA, nrow = 47, ncol = Nsim)
      sampCountyMod = matrix(NA, nrow = 47, ncol = Nsim)
      for(j in 1:Nsim){
        sampCounty[, j] = samp[[j]]$latent[1:47]
        if(includeCluster) {
          # if cluster effect is included, must debias predictions in each modeled strata
          clusterSigma = sqrt(1/samp[[j]]$hyperpar[3])
          muSigmaMat = cbind(sampCounty[, j], clusterSigma)
          sampCountyMod[, j] = logitNormMean(muSigmaMat = muSigmaMatUrban)
        }
      }
      sampCountyOverSampDat[,,i] = sampCounty
      sampCountyOverSampDatMod[,,i] = logit(sampCountyMod)
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
  Q10 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDat)[3])
  Q50 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDat)[3])
  Q90 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDat)[3])
  mm = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDat)[3])
  ss = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDat)[3])
  for(i in 1:dim(sampCountySRSDat)[3]){
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
  
  if(includeCluster) {
    Q10 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatMod)[3])
    Q50 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatMod)[3])
    Q90 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatMod)[3])
    mm = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatMod)[3])
    ss = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatMod)[3])
    for(i in 1:dim(sampCountySRSDatMod)[3]){
      tmp = processSamples(sampCountySRSDatMod[,,i])
      Q10[,i] = tmp$logit$CI[,1]
      Q50[,i] = tmp$logit$CI[,2]
      Q90[,i] = tmp$logit$CI[,3]
      mm[,i] = tmp$logit$mean
      ss[,i] = tmp$logit$stddev
    }
    resSRSdatMod = list(Q10 = Q10,
                     Q50 = Q50,
                     Q90 = Q90,
                     mean = mm,
                     stddev = ss)
  } else {
    resSRSdatMod = NULL
  }
  
  ## overSampDat
  Q10 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDat)[3])
  Q50 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDat)[3])
  Q90 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDat)[3])
  mm = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDat)[3])
  ss = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDat)[3])
  for(i in 1:dim(sampCountyOverSampDat)[3]){
    tmp = processSamples(sampCountyOverSampDat[,,i])
    Q10[,i] = tmp$logit$CI[,1]
    Q50[,i] = tmp$logit$CI[,2]
    Q90[,i] = tmp$logit$CI[,3]
    mm[,i] = tmp$logit$mean
    ss[,i] = tmp$logit$stddev
  }
  resOverSampDat = list(Q10 = Q10,
                       Q50 = Q50,
                       Q90 = Q90,
                       mean = mm,
                       stddev = ss)
  
  if(includeCluster) {
    Q10 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatMod)[3])
    Q50 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatMod)[3])
    Q90 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatMod)[3])
    mm = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatMod)[3])
    ss = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatMod)[3])
    for(i in 1:dim(sampCountyOverSampDatMod)[3]){
      tmp = processSamples(sampCountyOverSampDatMod[,,i])
      Q10[,i] = tmp$logit$CI[,1]
      Q50[,i] = tmp$logit$CI[,2]
      Q90[,i] = tmp$logit$CI[,3]
      mm[,i] = tmp$logit$mean
      ss[,i] = tmp$logit$stddev
    }
    resOverSampDatMod = list(Q10 = Q10,
                         Q50 = Q50,
                         Q90 = Q90,
                         mean = mm,
                         stddev = ss)
  } else {
    resOverSampDatMod = NULL
  }
  
  # Full result
  designRes = list(SRSdat = resSRSdat,
                   overSampDat = resOverSampDat)
  # save(file = 'kenyaSpatialDesignResultNew.RData', designRes = designRes)
  # save(file = paste0('kenyaSpatialDesignResultNewTausq0UrbRur', 
  #                      includeUrbanRural, '.RData'), designRes = designRes)
  testText = ifelse(test, "Test", "")
  save(file = paste0('kenyaSpatialDesignResultNewTausq', round(tausq, 4), 'UrbRur',
                     includeUrbanRural, 'Cluster', includeCluster, testText, '.RData'), 
       designRes = designRes)
  
  # include the debiased results if cluster effect is included
  if(includeCluster) {
    designRes = list(SRSdat = resSRSdat,
                     overSampDat = resOverSampDat)
    
    save(file = paste0('kenyaSpatialDesignResultNewTausq', round(tausq, 4), 'UrbRur',
                       includeUrbanRural, 'Cluster', includeCluster, 'debiased', testText, '.RData'), 
         designRes = designRes)
  }
  
  invisible(NULL)
}

