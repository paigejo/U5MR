##################################################################################
# designBased.R                                                                  #
#    Produce county level estimates using pixel-wise spatial model with          #
#    aggregation using true population and BYM spatial model                     #
##################################################################################

# same as runBYM, except fits the BYM2 reparameterized model (iid component in the BYM2 is that the county level)
runBYM2 = function(tausq=0.1^2, test=FALSE, includeUrbanRural=TRUE, includeCluster=TRUE, maxDataSets=NULL, 
                   aggregateByPopulation=FALSE, margVar=0.15^2, gamma=-1, plotPriorPost=FALSE) {
  
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
    parNames = c("Intercept", "Urban", "Cluster Var", "BYM2 Phi", "BYM2 Tot. Var", "BYM2 Spatial Var", "BYM2 iid Var")
  else
    parNames = c("Intercept", "Urban", "BYM2 Phi", "BYM2 Tot. Var", "BYM2 Spatial Var", "BYM2 iid Var")
  includeI = c(1, rep(2, includeUrbanRural), 3:length(parNames))
  parNames = parNames[includeI]
  
  # set the maximum number of data sets to analyze
  if(is.null(maxDataSets)) {
    maxDataSets = length(SRSDat$clustDat)
  }
  
  # Go through datasets for SRSDat
  sampCountySRSDat = array(NA, dim = c(47, Nsim, maxDataSets))
  sampCountySRSDatMod = array(NA, dim = c(47, Nsim, maxDataSets))
  
  # for the final parameters to store, 2 fixed effects, 2 + includeCluster estimated 
  # hyperparameters, and 2 hyperparameters we will get via transformation
  sampCountySRSDatPar = matrix(NA, 5 + includeUrbanRural + includeCluster, ncol=maxDataSets)
  sampCountySRSDatSD = matrix(NA, 5 + includeUrbanRural + includeCluster, ncol=maxDataSets)
  sampCountySRSDat10 = matrix(NA, 5 + includeUrbanRural + includeCluster, ncol=maxDataSets)
  sampCountySRSDat90 = matrix(NA, 5 + includeUrbanRural + includeCluster, ncol=maxDataSets)
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
                  control.compute = list(config = TRUE), 
                  quantiles=c(0.1, 0.5, 0.9))
    
    if(plotPriorPost) {
      # maxX = max(result$marginals.hyperpar[[1]][,1]) * 1.1
      maxX = inla.qmarginal(0.8, result$marginals.hyperpar[[1]])
      xs = exp(seq(-5, log(maxX), l=1500))
      ys = inla.pc.dprec(xs, 1, 0.01)
      maxY = max(max(ys), max(result$marginals.hyperpar[[1]][,2]))
      plot(xs, ys, type="l", col="blue", xlab="BYM2 Precision", ylab="Density", xlim=c(0, maxX), 
           ylim=c(0, maxY), main="BYM2 Precision Prior vs Posterior")
      lines(result$marginals.hyperpar[[1]], col="red")
      legend("topright", c("Prior", "Posterior"), col=c("blue", "red"), lty=1)
    }
    
    ## include parameter estimates in the table
    # fixed effects
    sampCountySRSDatPar[1:(1 + includeUrbanRural), i] = result$summary.fixed[,1]
    sampCountySRSDatSD[1:(1 + includeUrbanRural), i] = result$summary.fixed[,2]
    sampCountySRSDat10[1:(1 + includeUrbanRural), i] = result$summary.fixed[,3]
    sampCountySRSDat90[1:(1 + includeUrbanRural), i] = result$summary.fixed[,5]
    
    ## include parameter estimates in the table
    # fixed effects
    sampCountySRSDatPar[1:(1 + includeUrbanRural), i] = result$summary.fixed[,1]
    sampCountySRSDatSD[1:(1 + includeUrbanRural), i] = result$summary.fixed[,2]
    sampCountySRSDat10[1:(1 + includeUrbanRural), i] = result$summary.fixed[,3]
    sampCountySRSDat90[1:(1 + includeUrbanRural), i] = result$summary.fixed[,5]
    
    # BYM2 hyperparameter phi
    sampCountySRSDatPar[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,1]
    sampCountySRSDatSD[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,2]
    sampCountySRSDat10[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,3]
    sampCountySRSDat90[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,5]
    
    ## transformed hyperparameters
    # sample the hyperparameters, using the marginals to improve the sampling
    out = inla.hyperpar.sample(1000, result, improve.marginals=TRUE)
    transformFunction = function(x) {c(1/x[1], 1/x[1]*x[2], 1/x[1]*(1-x[2]), 1/x[3])}
    if(!includeCluster)
      transformFunction = function(x) {c(1/x[1], 1/x[1]*x[2], 1/x[1]*(1-x[2]))}
    transformedOut = apply(out, 1, transformFunction)
    
    # now calculate the summary statistics of the transformed BYM2 hyperparameters
    sampCountySRSDatPar[(3 + includeUrbanRural + includeCluster):(5 + includeUrbanRural + includeCluster), i] = rowMeans(transformedOut[1:3,])
    sampCountySRSDatSD[(3 + includeUrbanRural + includeCluster):(5 + includeUrbanRural + includeCluster), i] = apply(transformedOut[1:3,], 1, sd)
    sampCountySRSDat10[(3 + includeUrbanRural + includeCluster):(5 + includeUrbanRural + includeCluster), i] = apply(transformedOut[1:3,], 1, quantile, probs=.1)
    sampCountySRSDat90[(3 + includeUrbanRural + includeCluster):(5 + includeUrbanRural + includeCluster), i] = apply(transformedOut[1:3,], 1, quantile, probs=.9)
    
    # calculate summary statistics for cluster variance if necessary
    if(includeCluster) {
      sampCountySRSDatPar[2 + includeUrbanRural, i] = mean(transformedOut[4,])
      sampCountySRSDatSD[2 + includeUrbanRural, i] = sd(transformedOut[4,])
      sampCountySRSDat10[2 + includeUrbanRural, i] = quantile(transformedOut[4,], probs=.1)
      sampCountySRSDat90[2 + includeUrbanRural, i] = quantile(transformedOut[4,], probs=.9)
    }
    
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
  sampCountyOverSampDatPar = matrix(NA, 5 + includeUrbanRural + includeCluster, ncol=maxDataSets)
  sampCountyOverSampDatSD = matrix(NA, 5 + includeUrbanRural + includeCluster, ncol=maxDataSets)
  sampCountyOverSampDat10 = matrix(NA, 5 + includeUrbanRural + includeCluster, ncol=maxDataSets)
  sampCountyOverSampDat90 = matrix(NA, 5 + includeUrbanRural + includeCluster, ncol=maxDataSets)
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
    sampCountyOverSampDatPar[1:(1 + includeUrbanRural), i] = result$summary.fixed[,1]
    sampCountyOverSampDatSD[1:(1 + includeUrbanRural), i] = result$summary.fixed[,2]
    sampCountyOverSampDat10[1:(1 + includeUrbanRural), i] = result$summary.fixed[,3]
    sampCountyOverSampDat90[1:(1 + includeUrbanRural), i] = result$summary.fixed[,5]
    
    ## include parameter estimates in the table
    # fixed effects
    sampCountyOverSampDatPar[1:(1 + includeUrbanRural), i] = result$summary.fixed[,1]
    sampCountyOverSampDatSD[1:(1 + includeUrbanRural), i] = result$summary.fixed[,2]
    sampCountyOverSampDat10[1:(1 + includeUrbanRural), i] = result$summary.fixed[,3]
    sampCountyOverSampDat90[1:(1 + includeUrbanRural), i] = result$summary.fixed[,5]
    
    # BYM2 hyperparameter phi
    sampCountyOverSampDatPar[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,1]
    sampCountyOverSampDatSD[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,2]
    sampCountyOverSampDat10[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,3]
    sampCountyOverSampDat90[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,5]
    
    ## transformed hyperparameters
    # sample the hyperparameters, using the marginals to improve the sampling
    out = inla.hyperpar.sample(1000, result, improve.marginals=TRUE)
    transformFunction = function(x) {c(1/x[1], 1/x[1]*x[2], 1/x[1]*(1-x[2]), 1/x[3])}
    if(!includeCluster)
      transformFunction = function(x) {c(1/x[1], 1/x[1]*x[2], 1/x[1]*(1-x[2]))}
    transformedOut = apply(out, 1, transformFunction)
    
    # now calculate the summary statistics of the transformed BYM2 hyperparameters
    sampCountyOverSampDatPar[(3 + includeUrbanRural + includeCluster):(5 + includeUrbanRural + includeCluster), i] = rowMeans(transformedOut[1:3,])
    sampCountyOverSampDatSD[(3 + includeUrbanRural + includeCluster):(5 + includeUrbanRural + includeCluster), i] = apply(transformedOut[1:3,], 1, sd)
    sampCountyOverSampDat10[(3 + includeUrbanRural + includeCluster):(5 + includeUrbanRural + includeCluster), i] = apply(transformedOut[1:3,], 1, quantile, probs=.1)
    sampCountyOverSampDat90[(3 + includeUrbanRural + includeCluster):(5 + includeUrbanRural + includeCluster), i] = apply(transformedOut[1:3,], 1, quantile, probs=.9)
    
    # calculate summary statistics for cluster variance if necessary
    if(includeCluster) {
      sampCountyOverSampDatPar[2 + includeUrbanRural, i] = mean(transformedOut[4,])
      sampCountyOverSampDatSD[2 + includeUrbanRural, i] = sd(transformedOut[4,])
      sampCountyOverSampDat10[2 + includeUrbanRural, i] = quantile(transformedOut[4,], probs=.1)
      sampCountyOverSampDat90[2 + includeUrbanRural, i] = quantile(transformedOut[4,], probs=.9)
    }
    
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
  # make rural parameter the urban parameter
  if(includeUrbanRural) {
    sampCountySRSDatPar[2,] = -sampCountySRSDatPar[2,]
    sampCountySRSDat10[2,] = -sampCountySRSDat10[2,]
    sampCountySRSDat90[2,] = -sampCountySRSDat90[2,]
  }
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
  # make rural parameter the urban parameter
  if(includeUrbanRural) {
    sampCountyOverSampDatPar[2,] = -sampCountyOverSampDatPar[2,]
    sampCountyOverSampDat10[2,] = -sampCountyOverSampDat10[2,]
    sampCountyOverSampDat90[2,] = -sampCountyOverSampDat90[2,]
  }
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
                     includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 
                     maxDataSets, testText, '.RData'), 
       designRes = designRes)
  
  # include the debiased results if cluster effect is included
  if(includeCluster) {
    designRes = list(SRSdat = resSRSdatMod,
                     overSampDat = resOverSampDatMod, 
                     overSampDatPar = resOverSampDatPar, 
                     SRSdatPar = resSRSdatPar)
    
    save(file = paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                       includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, 'debiasedMaxDataSets', 
                       maxDataSets, testText, '.RData'), 
         designRes = designRes)
  }
  
  invisible(NULL)
}

# same as runBYM2, except fits a single data set (the ed global data frame)
runBYM2Dat = function(dat=ed, includeUrbanRural=TRUE, includeCluster=TRUE, saveResults=TRUE, fileNameRoot="Ed") {
  includeUrban = includeUrbanRural
  
  # Get true ratios of urban/rural
  urbRatio = vector('numeric', length = 47)
  counties = sort(unique(as.character(dat$admin1)))
  urbRatio = poppc$popUrb / poppc$popTotal
  sortI = matchMultiple(counties, poppc$County)
  urbRatio = urbRatio[sortI]
  
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
    parNames = c("Intercept", "Urban", "Cluster Var", "BYM2 Phi", "BYM2 Tot. Var", "BYM2 Spatial Var", "BYM2 iid Var")
  else
    parNames = c("Intercept", "Urban", "BYM2 Phi", "BYM2 Tot. Var", "BYM2 Spatial Var", "BYM2 iid Var")
  includeI = c(1, rep(2, includeUrban), 3:length(parNames))
  parNames = parNames[includeI]
  
  # Go through education data set
  sampCountyDat = matrix(NA, nrow=47, ncol=Nsim)
  
  # for the final parameters to store, 2 fixed effects, 2 + includeCluster estimated 
  # hyperparameters, and 2 hyperparameters we will get via transformation
  sampCountyDatPar = numeric(5 + includeUrban + includeCluster)
  sampCountyDatSD = numeric(5 + includeUrban + includeCluster)
  sampCountyDat10 = numeric(5 + includeUrban + includeCluster)
  sampCountyDat90 = numeric(5 + includeUrban + includeCluster)
  
  names(sampCountyDatPar) = parNames
  names(sampCountyDatSD) = parNames
  names(sampCountyDat10) = parNames
  names(sampCountyDat90) = parNames
  
  # Extract data
  dat$admin1 = factor(dat$admin1)
  
  # INLA data
  dat = list(y = dat$y,
             Ntrials = dat$n,
             rural = 1-dat$urban,
             idx = as.numeric(dat$admin1),
             idxEps = 1:length(dat$y))
  
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
                control.compute = list(config = TRUE), 
                quantiles=c(0.1, 0.5, 0.9))
  
  ## include parameter estimates in the table
  # fixed effects
  sampCountyDatPar[1:(1 + includeUrban)] = result$summary.fixed[,1]
  sampCountyDatSD[1:(1 + includeUrban)] = result$summary.fixed[,2]
  sampCountyDat10[1:(1 + includeUrban)] = result$summary.fixed[,3]
  sampCountyDat90[1:(1 + includeUrban)] = result$summary.fixed[,5]
  
  # BYM2 hyperparameter phi
  sampCountyDatPar[(2 + includeCluster + includeUrban)] = result$summary.hyperpar[2,1]
  sampCountyDatSD[(2 + includeCluster + includeUrban)] = result$summary.hyperpar[2,2]
  sampCountyDat10[(2 + includeCluster + includeUrban)] = result$summary.hyperpar[2,3]
  sampCountyDat90[(2 + includeCluster + includeUrban)] = result$summary.hyperpar[2,5]
  
  ## transformed hyperparameters
  # sample the hyperparameters, using the marginals to improve the sampling
  out = inla.hyperpar.sample(1000, result, improve.marginals=TRUE)
  transformFunction = function(x) {c(1/x[1], 1/x[1]*x[2], 1/x[1]*(1-x[2]), 1/x[3])}
  if(!includeCluster)
    transformFunction = function(x) {c(1/x[1], 1/x[1]*x[2], 1/x[1]*(1-x[2]))}
  transformedOut = apply(out, 1, transformFunction)
  
  # now calculate the summary statistics of the transformed BYM2 hyperparameters
  sampCountyDatPar[(3 + includeUrban + includeCluster):(5 + includeUrban + includeCluster)] = rowMeans(transformedOut[1:3,])
  sampCountyDatSD[(3 + includeUrban + includeCluster):(5 + includeUrban + includeCluster)] = apply(transformedOut[1:3,], 1, sd)
  sampCountyDat10[(3 + includeUrban + includeCluster):(5 + includeUrban + includeCluster)] = apply(transformedOut[1:3,], 1, quantile, probs=.1)
  sampCountyDat90[(3 + includeUrban + includeCluster):(5 + includeUrban + includeCluster)] = apply(transformedOut[1:3,], 1, quantile, probs=.9)
  
  # calculate summary statistics for cluster variance if necessary
  if(includeCluster) {
    sampCountyDatPar[2 + includeUrban] = mean(transformedOut[4,])
    sampCountyDatSD[2 + includeUrban] = sd(transformedOut[4,])
    sampCountyDat10[2 + includeUrban] = quantile(transformedOut[4,], probs=.1)
    sampCountyDat90[2 + includeUrban] = quantile(transformedOut[4,], probs=.9)
  }
  
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
    sampCountyDat = logit(expit(sampUrban)*urbRatio + expit(sampRural)*(1-urbRatio))
    sampCountyDatMod = logit(sampUrbanMod*urbRatio + sampRuralMod*(1-urbRatio))
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
    sampCountyDat = sampCounty
    sampCountyDatMod = logit(sampCountyMod)
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
  
  ## generate predictions
  Q10 = numeric(47)
  Q50 = numeric(47)
  Q90 = numeric(47)
  mm = numeric(47)
  ss = numeric(47)
  
  tmp = processSamples(sampCountyDat)
  Q10 = tmp$logit$CI[,1]
  Q50 = tmp$logit$CI[,2]
  Q90 = tmp$logit$CI[,3]
  mm = tmp$logit$mean
  ss = tmp$logit$stddev
  
  resDat = list(Q10 = Q10,
                 Q50 = Q50,
                 Q90 = Q90,
                 mean = mm,
                 stddev = ss)
  
  if(includeCluster) {
    tmp = processSamples(sampCountyDatMod)
    Q10 = tmp$logit$CI[,1]
    Q50 = tmp$logit$CI[,2]
    Q90 = tmp$logit$CI[,3]
    mm = tmp$logit$mean
    ss = tmp$logit$stddev
    resDatMod = list(Q10 = Q10,
                      Q50 = Q50,
                      Q90 = Q90,
                      mean = mm,
                      stddev = ss)
  } else {
    resDatMod = NULL
  }
  
  ## now collect the parameters
  # make rural parameter the urban parameter
  if(includeUrban) {
    sampCountyDatPar[2] = -sampCountyDatPar[2]
    sampCountyDat10[2] = -sampCountyDat10[2]
    sampCountyDat90[2] = -sampCountyDat90[2]
  }
  mm = sampCountyDatPar
  ss = sampCountyDatSD
  Q10 = sampCountyDat10
  Q90 = sampCountyDat90
  resDatPar = data.frame(list(Q10 = Q10,
                               Q90 = Q90,
                               mean = mm,
                               stddev = ss))
  
  # Full result
  designRes = list(predictions = resDat,
                   parameters = resDatPar)
  # save(file = 'kenyaSpatialDesignResultNew.RData', designRes = designRes)
  # save(file = paste0('kenyaSpatialDesignResultNewTausq0UrbRur', 
  #                      includeUrbanRural, '.RData'), designRes = designRes)
  if(saveResults) {
    save(file = paste0('bym2', fileNameRoot, 'UrbRur',includeUrbanRural, 'Cluster', includeCluster, '.RData'), 
         designRes = designRes)
  }
  
  # include the debiased results if cluster effect is included
  if(includeCluster) {
    designRes = list(predictions = resDatMod,
                     parameters = resDatPar)
    
    if(saveResults) {
      save(file = paste0('bym2', fileNameRoot, 'UrbRur',includeUrbanRural, 'Cluster', includeCluster, 'debiased.RData'), 
           designRes = designRes)
    }
  }
  
  designRes = list(predictions = resDatMod,
                   parameters = resDatPar)
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

