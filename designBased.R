##################################################################################
# designBased.R                                                                  #
#    Produce county level estimates using pixel-wise spatial model with          #
#    aggregation using true population and BYM spatial model                     #
##################################################################################

# same as runBYM, except fits the BYM2 reparameterized model (iid component in the BYM2 is that the county level)
runBYM2 = function(tausq=0.1^2, test=FALSE, includeUrbanRural=TRUE, includeCluster=TRUE, maxDataSets=NULL, 
                   aggregateByPopulation=FALSE, margVar=0.15^2, gamma=-1) {
  
  # load and relevant data
  if(!test)
    load(paste0("simDataMultiBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                "HHoldVar0urbanOverSamplefrac0.RData"))
  else
    load(paste0("simDataMultiBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                "HHoldVar0urbanOverSamplefrac0Test.RData"))
  
  # Get true ratios of urban/rural
  urbRatio = vector('numeric', length = 47)
  counties = sort(unique(as.character(SRSDat$eaDat$admin1)))
  if(!aggregateByPopulation) {
    for(i in 1:47){
      idx = which(as.numeric(factor(SRSDat$eaDat$admin1)) == i)
      urbanI = idx[SRSDat$eaDat$urban[idx]]
      urbRatio[i] = sum(SRSDat$eaDat$numChildren[urbanI])/sum(SRSDat$eaDat$numChildren[idx])
    }
  } else {
    urbRatio = poppc$popUrb / poppc$popTotal
    sortI = matchMultiple(counties, poppc$County)
    urbRatio = urbRatio[sortI]
  }
  
  # Define formula
  if(includeUrbanRural) {
    if(includeCluster) {
      formula = y ~ rural +
        f(idx, model="bym2",
          graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
          hyper=list(prec=list(param=c(1, 0.01), prior="pc.prec"), phi=list(param=c(0.5, 0.5), prior="pc"))) +
        f(idxEps, model = "iid",
          hyper = list(prec = list(prior = "pc.prec", param = c(3,0.01))))
    } else {
      formula = y ~ rural + 
        f(idx, model="bym2",
          graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
          hyper=list(prec=list(param=c(1, 0.01), prior="pc.prec"), phi=list(param=c(0.5, 0.5), prior="pc")))
    }
  } else {
    if(includeCluster) {
      formula = y ~ f(idx, model="bym2",
                      graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
                      hyper=list(prec=list(param=c(1, 0.01), prior="pc.prec"), 
                                 phi=list(param=c(0.5, 0.5), prior="pc"))) +
        f(idxEps, model = "iid",
          hyper = list(prec = list(prior = "pc.prec", param = c(3,0.01))))
    } else {
      formula = y ~ 
        f(idx, model="bym2",
          graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
          hyper=list(prec=list(param=c(1, 0.01), prior="pc.prec"), phi=list(param=c(0.5, 0.5), prior="pc")))
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
  
  if(includeCluster)
    parNames = c("Intercept", "Urban", "BYM2 Tot. Var", "BYM2 Phi", "Cluster Var", "BYM2 Spatial Var", "BYM2 iid Var")
  else
    parNames = c("Intercept", "Urban", "BYM2 Tot. Var", "BYM2 Phi", "BYM2 Spatial Var", "BYM2 iid Var")
  
  # set the maximum number of data sets to analyze
  if(is.null(maxDataSets)) {
    maxDataSets = length(SRSDat$clustDat)
  }
  
  # Go through datasets for SRSDat
  sampCountySRSDat = array(NA, dim = c(47, Nsim, maxDataSets))
  sampCountySRSDatMod = array(NA, dim = c(47, Nsim, maxDataSets))
  
  # for the final parameters to store, 2 fixed effects, 2 + includeCluster estimated 
  # hyperparameters, and 2 hyperparameters we will get via transformation
  sampCountySRSDatPar = matrix(NA, nrow = 2 + 2 + includeCluster + 2, ncol=maxDataSets)
  sampCountySRSDatSD = matrix(NA, nrow = 2 + 2 + includeCluster + 2, ncol=maxDataSets)
  sampCountySRSDat10 = matrix(NA, nrow = 2 + 2 + includeCluster + 2, ncol=maxDataSets)
  sampCountySRSDat90 = matrix(NA, nrow = 2 + 2 + includeCluster + 2, ncol=maxDataSets)
  rownames(sampCountySRSDatPar) = parNames
  rownames(sampCountySRSDatSD) = parNames
  rownames(sampCountySRSDat10) = parNames
  rownames(sampCountySRSDat90) = parNames
  for(i in 1:maxDataSets){
    # Extract data
    currData = SRSDat$clustDat[[i]]
    currData$admin1 = factor(currData$admin1)
    
    # INLA data
    dat = list(y = currData$died,
               Ntrials = currData$numChildren,
               rural = 1-currData$urban,
               idx = as.numeric(currData$admin1),
               idxEps = 1:length(currData$died))
    
    # Add unobserved data to make sampling easier
    dat$y = c(rep(NA, 47*2), dat$y)
    dat$Ntrials = c(rep(1, 47*2), dat$Ntrials)
    dat$rural = c(rep(c(0,1), each = 47), dat$rural)
    dat$idx = c(rep(1:47, 2), dat$idx)
    dat$idxEps = c(rep(NA, 47*2), dat$idxEps)
    
    # Run model
    print(paste0("fitting SRS model for dataset ", i, "/", maxDataSets))
    result = inla(formula = formula, 
                  family="binomial",
                  Ntrials = Ntrials,
                  data=dat, 
                  control.compute = list(config = TRUE))
    
    ## include parameter estimates in the table
    # fixed effects
    sampCountySRSDatPar[1:2, i] = result$summary.fixed[,1]
    sampCountySRSDatSD[1:2, i] = result$summary.fixed[,2]
    sampCountySRSDat10[1:2, i] = result$summary.fixed[,3]
    sampCountySRSDat90[1:2, i] = result$summary.fixed[,5]
    
    # hyperparameters
    sampCountySRSDatPar[3:(4 + includeCluster), i] = result$summary.hyperpar[,1]
    sampCountySRSDatSD[3:(4 + includeCluster), i] = result$summary.hyperpar[,2]
    sampCountySRSDat10[3:(4 + includeCluster), i] = result$summary.hyperpar[,3]
    sampCountySRSDat90[3:(4 + includeCluster), i] = result$summary.hyperpar[,5]
    
    ## transformed hyperparameters
    # sample the hyperparameters, using the marginals to improve the sampling
    out = inla.hyperpar.sample(1000, result, improve.marginals=TRUE)
    transformedOut = apply(out, 1, function(x) {c(1/x[1]*x[2], 1/x[1]*(1-x[2]))})
    
    # now calculate their summary statistics
    sampCountySRSDatPar[(5 + includeCluster):(6 + includeCluster), i] = rowMeans(transformedOut)
    sampCountySRSDatSD[(5 + includeCluster):(6 + includeCluster), i] = apply(transformedOut, 1, sd)
    sampCountySRSDat10[(5 + includeCluster):(6 + includeCluster), i] = apply(transformedOut, 1, quantile, probs=.1)
    sampCountySRSDat90[(5 + includeCluster):(6 + includeCluster), i] = apply(transformedOut, 1, quantile, probs=.9)
    
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
  sampCountyOverSampDat = array(NA, dim = c(47, Nsim, maxDataSets))
  sampCountyOverSampDatMod = array(NA, dim = c(47, Nsim, maxDataSets))
  
  # for the final parameters to store, 2 fixed effects, 2 + includeCluster estimated 
  # hyperparameters, and 2 hyperparameters we will get via transformation
  sampCountyOverSampDatPar = matrix(NA, nrow = 2 + 2 + includeCluster + 2, ncol=maxDataSets)
  sampCountyOverSampDatSD = matrix(NA, nrow = 2 + 2 + includeCluster + 2, ncol=maxDataSets)
  sampCountyOverSampDat10 = matrix(NA, nrow = 2 + 2 + includeCluster + 2, ncol=maxDataSets)
  sampCountyOverSampDat90 = matrix(NA, nrow = 2 + 2 + includeCluster + 2, ncol=maxDataSets)
  rownames(sampCountyOverSampDatPar) = parNames
  rownames(sampCountyOverSampDatSD) = parNames
  rownames(sampCountyOverSampDat10) = parNames
  rownames(sampCountyOverSampDat90) = parNames
  for(i in 1:maxDataSets) {
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
    print(paste0("fitting urban oversampled model for dataset ", i, "/", maxDataSets))
    result = inla(formula = formula,
                  family="binomial",
                  Ntrials = Ntrials,
                  data=dat, quantiles = c(0.1, 0.5, 0.9), 
                  control.compute = list(config = TRUE))
    
    # fixed effects
    sampCountyOverSampDatPar[1:2, i] = result$summary.fixed[,1]
    sampCountyOverSampDatSD[1:2, i] = result$summary.fixed[,2]
    sampCountyOverSampDat10[1:2, i] = result$summary.fixed[,3]
    sampCountyOverSampDat90[1:2, i] = result$summary.fixed[,5]
    
    # hyperparameters
    sampCountyOverSampDatPar[3:(4 + includeCluster), i] = result$summary.hyperpar[,1]
    sampCountyOverSampDatSD[3:(4 + includeCluster), i] = result$summary.hyperpar[,2]
    sampCountyOverSampDat10[3:(4 + includeCluster), i] = result$summary.hyperpar[,3]
    sampCountyOverSampDat90[3:(4 + includeCluster), i] = result$summary.hyperpar[,5]
    
    ## transformed hyperparameters
    # sample the hyperparameters, using the marginals to improve the sampling
    out = inla.hyperpar.sample(1000, result, improve.marginals=TRUE)
    transformedOut = apply(out, 1, function(x) {c(1/x[1]*x[2], 1/x[1]*(1-x[2]))})
    
    # now calculate their summary statistics
    sampCountyOverSampDatPar[(5 + includeCluster):(6 + includeCluster), i] = rowMeans(transformedOut)
    sampCountyOverSampDatSD[(5 + includeCluster):(6 + includeCluster), i] = apply(transformedOut, 1, sd)
    sampCountyOverSampDat10[(5 + includeCluster):(6 + includeCluster), i] = apply(transformedOut, 1, quantile, probs=.1)
    sampCountyOverSampDat90[(5 + includeCluster):(6 + includeCluster), i] = apply(transformedOut, 1, quantile, probs=.9)
    
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
          sampCountyMod[, j] = logitNormMean(muSigmaMat = muSigmaMat)
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
  
  ## now collect the parameters
  # first invert precision to variance, make rural parameter the urban parameter
  precIds = 3
  if(includeCluster)
    precIds = c(precIds, 5)
  sampCountySRSDatPar[precIds, ] = 1 / sampCountySRSDatPar[precIds, ]
  sampCountySRSDatSD[precIds, ] = 1 / sampCountySRSDatSD[precIds, ]
  sampCountySRSDat10[precIds, ] = 1 / sampCountySRSDat10[precIds, ]
  sampCountySRSDat90[precIds, ] = 1 / sampCountySRSDat90[precIds, ]
  sampCountySRSDatPar[2, ] = -sampCountySRSDatPar[2, ]
  sampCountySRSDatSD[2, ] = -sampCountySRSDatSD[2, ]
  sampCountySRSDat10[2, ] = -sampCountySRSDat10[2, ]
  sampCountySRSDat90[2, ] = -sampCountySRSDat90[2, ]
  mm = rowMeans(sampCountySRSDatPar)
  ss = rowMeans(sampCountySRSDatSD)
  Q10 = rowMeans(sampCountySRSDat10)
  Q90 = rowMeans(sampCountySRSDat90)
  resSRSdatPar = data.frame(list(Q10 = Q10,
                                 Q90 = Q90,
                                 mean = mm,
                                 stddev = ss))
  
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
  
  ## now collect the parameters
  # first invert precision to variance, make rural parameter the urban parameter
  precIds = 3
  if(includeCluster)
    precIds = c(precIds, 5)
  sampCountyOverSampDatPar[precIds, ] = 1 / sampCountyOverSampDatPar[precIds, ]
  sampCountyOverSampDatSD[precIds, ] = 1 / sampCountyOverSampDatSD[precIds, ]
  sampCountyOverSampDat10[precIds, ] = 1 / sampCountyOverSampDat10[precIds, ]
  sampCountyOverSampDat90[precIds, ] = 1 / sampCountyOverSampDat90[precIds, ]
  sampCountyOverSampDatPar[2, ] = -sampCountyOverSampDatPar[2, ]
  sampCountyOverSampDatSD[2, ] = -sampCountyOverSampDatSD[2, ]
  sampCountyOverSampDat10[2, ] = -sampCountyOverSampDat10[2, ]
  sampCountyOverSampDat90[2, ] = -sampCountyOverSampDat90[2, ]
  mm = rowMeans(sampCountyOverSampDatPar)
  ss = rowMeans(sampCountyOverSampDatSD)
  Q10 = rowMeans(sampCountyOverSampDat10)
  Q90 = rowMeans(sampCountyOverSampDat90)
  resOverSampDatPar = data.frame(list(Q10 = Q10,
                                      Q90 = Q90,
                                      mean = mm,
                                      stddev = ss))
  
  # Full result
  designRes = list(SRSdat = resSRSdat,
                   overSampDat = resOverSampDat, 
                   overSampDatPar = resOverSampDatPar, 
                   SRSdatPar = resSRSdatPar)
  # save(file = 'kenyaSpatialDesignResultNew.RData', designRes = designRes)
  # save(file = paste0('kenyaSpatialDesignResultNewTausq0UrbRur', 
  #                      includeUrbanRural, '.RData'), designRes = designRes)
  
  testText = ifelse(test, "Test", "")
  save(file = paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                     includeUrbanRural, 'Cluster', includeCluster, "maxDataSets", maxDataSets, testText, '.RData'), 
       designRes = designRes)
  
  # include the debiased results if cluster effect is included
  if(includeCluster) {
    designRes = list(SRSdat = resSRSdatMod,
                     overSampDat = resOverSampDatMod, 
                     overSampDatPar = resOverSampDatPar, 
                     SRSdatPar = resSRSdatPar)
    
    save(file = paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                       includeUrbanRural, 'Cluster', includeCluster, 'debiasedMaxDataSets', maxDataSets, testText, '.RData'), 
         designRes = designRes)
  }
  
  invisible(NULL)
}

# same as runBYM2, except fits a single data set (the mort global data frame)
runBYM2Mort = function(dat=mort, includeUrbanRural=TRUE, includeCluster=TRUE) {
  
  # Get true ratios of urban/rural
  urbRatio = vector('numeric', length = 47)
  for(i in 1:47){
    # TODO: fix this!!!!!!!!!!!!!!!!!!!!
    idx = dat$admin1 == i
    urbanI = idx[SRSDat$eaDat$urban[idx]]
    urbRatio[i] = sum(SRSDat$eaDat$numChildren[urbanI])/sum(SRSDat$eaDat$numChildren[idx])
  }
  
  # Define formula
  if(includeUrbanRural) {
    if(includeCluster) {
      formula = y ~ rural +
        f(idx, model="bym2",
          graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
          hyper=list(prec=list(param=c(1, 0.01), prior="pc.prec"), phi=list(param=c(0.5, 0.5), prior="pc"))) +
        f(idxEps, model = "iid",
          hyper = list(prec = list(prior = "pc.prec", param = c(3,0.01))))
    } else {
      formula = y ~ rural + 
        f(idx, model="bym2",
          graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
          hyper=list(prec=list(param=c(1, 0.01), prior="pc.prec"), phi=list(param=c(0.5, 0.5), prior="pc")))
    }
  } else {
    if(includeCluster) {
      formula = y ~ f(idx, model="bym2",
                      graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
                      hyper=list(prec=list(param=c(1, 0.01), prior="pc.prec"), 
                                 phi=list(param=c(0.5, 0.5), prior="pc"))) +
        f(idxEps, model = "iid",
          hyper = list(prec = list(prior = "pc.prec", param = c(3,0.01))))
    } else {
      formula = y ~ 
        f(idx, model="bym2",
          graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
          hyper=list(prec=list(param=c(1, 0.01), prior="pc.prec"), phi=list(param=c(0.5, 0.5), prior="pc")))
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
  
  if(includeCluster)
    parNames = c("Intercept", "Urban", "BYM2 Tot. Var", "BYM2 Phi", "Cluster Var", "BYM2 Spatial Var", "BYM2 iid Var")
  else
    parNames = c("Intercept", "Urban", "BYM2 Tot. Var", "BYM2 Phi", "BYM2 Spatial Var", "BYM2 iid Var")
  
  # Go through datasets for SRSDat
  sampCountySRSDat = matrix(NA, nrow=47, ncol=Nsim)
  sampCountySRSDatMod = matrix(NA, nrow=47, ncol=Nsim)
  
  # for the final parameters to store, 2 fixed effects, 2 + includeCluster estimated 
  # hyperparameters, and 2 hyperparameters we will get via transformation
  sampCountySRSDatPar = numeric(2 + 2 + includeCluster + 2)
  sampCountySRSDatSD = numeric(2 + 2 + includeCluster + 2)
  sampCountySRSDat10 = numeric(2 + 2 + includeCluster + 2)
  sampCountySRSDat90 = numeric(2 + 2 + includeCluster + 2)
  names(sampCountySRSDatPar) = parNames
  names(sampCountySRSDatSD) = parNames
  names(sampCountySRSDat10) = parNames
  names(sampCountySRSDat90) = parNames
  
    # Extract data
    mort$admin1 = factor(mort$admin1)
    
    # INLA data
    dat = list(y = mort$died,
               Ntrials = mort$numChildren,
               rural = 1-mort$urban,
               idx = as.numeric(mort$admin1),
               idxEps = 1:length(mort$died))
    
    # Add unobserved data to make sampling easier
    dat$y = c(rep(NA, 47*2), dat$y)
    dat$Ntrials = c(rep(1, 47*2), dat$Ntrials)
    dat$rural = c(rep(c(0,1), each = 47), dat$rural)
    dat$idx = c(rep(1:47, 2), dat$idx)
    dat$idxEps = c(rep(NA, 47*2), dat$idxEps)
    
    # Run model
    print("fitting BYM model...")
    result = inla(formula = formula, 
                  family="binomial",
                  Ntrials = Ntrials,
                  data=dat, 
                  control.compute = list(config = TRUE))
    
    ## include parameter estimates in the table
    # fixed effects
    sampCountySRSDatPar[1:2] = result$summary.fixed[,1]
    sampCountySRSDatSD[1:2] = result$summary.fixed[,2]
    sampCountySRSDat10[1:2] = result$summary.fixed[,3]
    sampCountySRSDat90[1:2] = result$summary.fixed[,5]
    
    # hyperparameters
    sampCountySRSDatPar[3:(4 + includeCluster)] = result$summary.hyperpar[,1]
    sampCountySRSDatSD[3:(4 + includeCluster)] = result$summary.hyperpar[,2]
    sampCountySRSDat10[3:(4 + includeCluster)] = result$summary.hyperpar[,3]
    sampCountySRSDat90[3:(4 + includeCluster)] = result$summary.hyperpar[,5]
    
    ## transformed hyperparameters
    # sample the hyperparameters, using the marginals to improve the sampling
    out = inla.hyperpar.sample(1000, result, improve.marginals=TRUE)
    transformedOut = apply(out, 1, function(x) {c(1/x[1]*x[2], 1/x[1]*(1-x[2]))})
    
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
      sampCountySRSDat = logit(expit(sampUrban)*urbRatio + expit(sampRural)*(1-urbRatio))
      sampCountySRSDatMod = logit(sampUrbanMod*urbRatio + sampRuralMod*(1-urbRatio))
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
      sampCountySRSDat = sampCounty
      sampCountySRSDatMod = logit(sampCountyMod)
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
  Q10 = numeric(47)
  Q50 = numeric(47)
  Q90 = numeric(47)
  mm = numeric(47)
  ss = numeric(47)
  
  tmp = processSamples(sampCountySRSDat)
  Q10 = tmp$logit$CI[,1]
  Q50 = tmp$logit$CI[,2]
  Q90 = tmp$logit$CI[,3]
  mm = tmp$logit$mean
  ss = tmp$logit$stddev
    
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
  
  ## now collect the parameters
  # first invert precision to variance, make rural parameter the urban parameter
  precIds = 3
  if(includeCluster)
    precIds = c(precIds, 5)
  sampCountySRSDatPar[precIds, ] = 1 / sampCountySRSDatPar[precIds, ]
  sampCountySRSDatSD[precIds, ] = 1 / sampCountySRSDatSD[precIds, ]
  sampCountySRSDat10[precIds, ] = 1 / sampCountySRSDat10[precIds, ]
  sampCountySRSDat90[precIds, ] = 1 / sampCountySRSDat90[precIds, ]
  sampCountySRSDatPar[2, ] = -sampCountySRSDatPar[2, ]
  sampCountySRSDatSD[2, ] = -sampCountySRSDatSD[2, ]
  sampCountySRSDat10[2, ] = -sampCountySRSDat10[2, ]
  sampCountySRSDat90[2, ] = -sampCountySRSDat90[2, ]
  mm = rowMeans(sampCountySRSDatPar)
  ss = rowMeans(sampCountySRSDatSD)
  Q10 = rowMeans(sampCountySRSDat10)
  Q90 = rowMeans(sampCountySRSDat90)
  resSRSdatPar = data.frame(list(Q10 = Q10,
                                 Q90 = Q90,
                                 mean = mm,
                                 stddev = ss))
  
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
  
  ## now collect the parameters
  # first invert precision to variance, make rural parameter the urban parameter
  precIds = 3
  if(includeCluster)
    precIds = c(precIds, 5)
  sampCountyOverSampDatPar[precIds, ] = 1 / sampCountyOverSampDatPar[precIds, ]
  sampCountyOverSampDatSD[precIds, ] = 1 / sampCountyOverSampDatSD[precIds, ]
  sampCountyOverSampDat10[precIds, ] = 1 / sampCountyOverSampDat10[precIds, ]
  sampCountyOverSampDat90[precIds, ] = 1 / sampCountyOverSampDat90[precIds, ]
  sampCountyOverSampDatPar[2, ] = -sampCountyOverSampDatPar[2, ]
  sampCountyOverSampDatSD[2, ] = -sampCountyOverSampDatSD[2, ]
  sampCountyOverSampDat10[2, ] = -sampCountyOverSampDat10[2, ]
  sampCountyOverSampDat90[2, ] = -sampCountyOverSampDat90[2, ]
  mm = rowMeans(sampCountyOverSampDatPar)
  ss = rowMeans(sampCountyOverSampDatSD)
  Q10 = rowMeans(sampCountyOverSampDat10)
  Q90 = rowMeans(sampCountyOverSampDat90)
  resOverSampDatPar = data.frame(list(Q10 = Q10,
                                      Q90 = Q90,
                                      mean = mm,
                                      stddev = ss))
  
  # Full result
  designRes = list(SRSdat = resSRSdat,
                   overSampDat = resOverSampDat, 
                   overSampDatPar = resOverSampDatPar, 
                   SRSdatPar = resSRSdatPar)
  # save(file = 'kenyaSpatialDesignResultNew.RData', designRes = designRes)
  # save(file = paste0('kenyaSpatialDesignResultNewTausq0UrbRur', 
  #                      includeUrbanRural, '.RData'), designRes = designRes)
  
  testText = ifelse(test, "Test", "")
  save(file = paste0('bym2Tausq', round(tausq, 4), 'UrbRur',
                     includeUrbanRural, 'Cluster', includeCluster, "maxDataSets", maxDataSets, testText, '.RData'), 
       designRes = designRes)
  
  # include the debiased results if cluster effect is included
  if(includeCluster) {
    designRes = list(SRSdat = resSRSdatMod,
                     overSampDat = resOverSampDatMod, 
                     overSampDatPar = resOverSampDatPar, 
                     SRSdatPar = resSRSdatPar)
    
    save(file = paste0('bym2Tausq', round(tausq, 4), 'UrbRur',
                       includeUrbanRural, 'Cluster', includeCluster, 'debiasedMaxDataSets', maxDataSets, testText, '.RData'), 
         designRes = designRes)
  }
  
  invisible(NULL)
}

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
  sampCountySRSDat = array(NA, dim = c(47, Nsim, maxDataSets))
  sampCountySRSDatMod = array(NA, dim = c(47, Nsim, maxDataSets))
  for(i in 1:maxDataSets){
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
  sampCountyOverSampDat = array(NA, dim = c(47, Nsim, maxDataSets))
  sampCountyOverSampDatMod = array(NA, dim = c(47, Nsim, maxDataSets))
  for(i in 1:maxDataSets) {
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
          sampCountyMod[, j] = logitNormMean(muSigmaMat = muSigmaMat)
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
    designRes = list(SRSdat = resSRSdatMod,
                     overSampDat = resOverSampDatMod)
    
    save(file = paste0('kenyaSpatialDesignResultNewTausq', round(tausq, 4), 'UrbRur',
                       includeUrbanRural, 'Cluster', includeCluster, 'debiased', testText, '.RData'), 
         designRes = designRes)
  }
  
  invisible(NULL)
}

