# code for generating table for Jon's talk with model scores

# generate predictive results given the simulated datasets simDataMulti.RData
# includeClustEffect, urbanEffect: control whether or knot to include a 
# cluster random effect in urban fixed effect in the model
# genRegionLevel, keepPixelPreds, genEALevel: control whether to include 
# regional, pixel, and enumeration area level predictions
resultsSPDE = function(nPostSamples=100, test=FALSE, nTest=5, verbose=TRUE, 
                       includeClustEffect=FALSE, int.strategy="eb", 
                       genRegionLevel=TRUE, keepPixelPreds=TRUE, 
                       genEALevel=TRUE, urbanEffect=TRUE, tausq=0, 
                       genCountLevel=FALSE, saveResults=!test, exactAggregation=genCountLevel) {
  # Load data
  # load("simDataMulti.RData") # overSampDat, SRSDat
  # load a different 1 of these depending on whether a cluster effect should be included 
  # in the simulation of the data or not (tausq is the cluster effect variance)
  # load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOver2.RData")
  # load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOver2.RData")
  load(paste0("simDataMultiBeta-1.75margVar0.0225tausq", round(tausq, 3), 
              "gamma-1HHoldVar0urbanOver2.RData"))
  eaDat = overSampDat$eaDat
  clustSRS = SRSDat$clustDat
  clustOverSamp = overSampDat$clustDat
  
  if(test) {
    clustSRS = lapply(1:nTest, function(i) {clustSRS[[i]]})
    clustOverSamp = lapply(1:nTest, function(i) {clustOverSamp[[i]]})
  }
  
  # SRS results are correct, so load those and recompute overSamp results
  # load(paste0("resultsSPDETausq", round(tausq, 4), 
  #             "urbanEffect", as.character(urbanEffect), ".RData"))
  spdeSRS = resultsSPDEHelper(clustSRS, eaDat, nPostSamples = nPostSamples, verbose=verbose,
                              includeClustEffect=includeClustEffect, int.strategy=int.strategy,
                              genRegionLevel=genRegionLevel, keepPixelPreds=keepPixelPreds,
                              genEALevel=genEALevel, urbanEffect=urbanEffect, 
                              genCountLevel=genCountLevel, exactAggregation=exactAggregation)
  spdeOverSamp = resultsSPDEHelper(clustOverSamp, eaDat, nPostSamples = nPostSamples, verbose=verbose, 
                                   includeClustEffect=includeClustEffect, int.strategy=int.strategy, 
                                   genRegionLevel=genRegionLevel, keepPixelPreds=keepPixelPreds, 
                                   genEALevel=genEALevel, urbanEffect=urbanEffect, 
                                   genCountLevel=genCountLevel, exactAggregation=exactAggregation)
  
  if(saveResults) {
    save(spdeSRS, spdeOverSamp, file=paste0("resultsSPDETausq", round(tausq, 4), 
                                            "urbanEffect", as.character(urbanEffect), ".RData"))
  }
  
  
  list(spdeSRS=spdeSRS, spdeOverSamp=spdeOverSamp)
}

# generate predictive results given a simulated dataset simDataMulti.RData
resultsSPDEHelper = function(clustDatMulti, eaDat, nPostSamples=100, verbose=TRUE, 
                             includeClustEffect=FALSE, int.strategy="eb", 
                             genRegionLevel=TRUE, keepPixelPreds=TRUE, genEALevel=TRUE, 
                             urbanEffect=TRUE, genCountLevel=FALSE, exactAggregation=genCountLevel, 
                             predictionType=c("mean", "median"), parClust=cl) {
  
  # match the requested prediction type with one of the possible options
  predictionType = match.arg(predictionType)
  
  # get the true population mortality rate for each county
  out = getTruthByCounty(eaDat)
  counties = out$counties
  trueMort = out$mortRate
  nsim = length(clustDatMulti)
  
  # get truth
  counties = as.character(unique(mort$admin1))
  regions=as.character(unique(mort$region))
  truth = getTruthByCounty(eaDat, counties)
  trueMort = truth$mortRate
  
  # get prediction locations from population grid
  # popGrid = makeInterpPopGrid()
  load("popGrid.RData")
  predCoords = cbind(popGrid$east, popGrid$north)
  predUrban = popGrid$urban
  if(genEALevel) {
    # Must predict at enumeration areas as well. Include enumeration areas as 
    # first rows of prediction coordinates an prediction urban/rural
    predCoords = rbind(cbind(eaDat$east, eaDat$north), predCoords)
    predUrban = c(eaDat$urban, predUrban)
  }
  # we only care about the probability, not counts, so not used except for the purposes 
  # of calling inla:
  predNs = rep(25, nrow(predCoords))
  
  # compute Bias & MSE & mean(Var) & 80\% coverage for each simulation
  if(is.null(parClust)) {
    countyResults = as.list(1:nsim)
    regionResults = as.list(1:nsim)
    pixelResults = as.list(1:nsim)
    eaResults = as.list(1:nsim)
    for(i in 1:nsim) {
      print(paste0("iteration ", i, "/", nsim))
      
      # get the simulated sample
      clustDat = clustDatMulti[[i]]
      
      # get observations from dataset
      obsCoords = cbind(clustDat$east, clustDat$north)
      obsNs = rep(25, nrow(obsCoords))
      obsCounts = clustDat$died
      obsUrban = clustDat$urban
      
      # fit model, get all predictions for each areal level and each posterior sample
      fit = fitSPDEModel(obsCoords, obsNs=obsNs, obsCounts, obsUrban, predCoords, predNs=predNs, 
                         predUrban, genCountyLevel=TRUE, popGrid=popGrid, nPostSamples=nPostSamples, 
                         verbose = verbose, includeClustEffect=includeClustEffect, 
                         int.strategy=int.strategy, genRegionLevel=genRegionLevel, 
                         keepPixelPreds=keepPixelPreds, genEALevel=genEALevel, 
                         urbanEffect=urbanEffect, link=1, predictionType=predictionType, 
                         exactAggregation=exactAggregation, genCountLevel=genCountLevel, 
                         eaDat=eaDat)
      countyPredMat = fit$countyPredMat
      regionPredMat = fit$regionPredMat
      pixelPredMat = fit$pixelPredMat
      eaPredMat = fit$eaPredMat
      eaMarginals = fit$eaMarginals
      
      ## calculate model fit properties for each level of predictions
      
      # enumeration area level
      if(genEALevel) {
        if(!genCountLevel) {
          if(predictionType == "mean")
            thisu1mEA = rowMeans(eaPredMat)
          else
            thisu1mEA = apply(eaPredMat, 1, median)
          thislowerEA = logit(apply(eaPredMat, 1, function(x) {quantile(x, probs=.1)}))
          thisupperEA = logit(apply(eaPredMat, 1, function(x) {quantile(x, probs=.9)}))
          thisvar.estEA = apply(logit(eaPredMat), 1, var)
          thisRightRejectEA = NA
          thisLeftRejectEA = NA
        }
        else {
          # get the marginal "binomial" densities at each location
          allMarginals = fit$mod$marginals.linear.predictor[fit$predInds]
          n = 25
          binomProb = function(k, i) {
            inla.emarginal(function(logitP) {dbinom(k, n, expit(logitP))}, allMarginals[[i]])
          }
          
          ## make highest density coverage interval on count scale
          # generate the "binomial" pmfs for each marginal
          probs = matrix(mapply(binomProb, 
                                k = matrix(0:n, nrow=length(allMarginals), ncol=n + 1, byrow = TRUE), 
                                i = matrix(1:length(allMarginals), nrow=length(allMarginals), ncol=n + 1)), 
                         nrow=length(allMarginals), ncol=n + 1)
          
          # now generate the summary statistics about the marginals
          if(predictionType == "mean")
            thisu1mEA = rowSums(sweep(probs, 2, STATS = 0:n, "*"))
          else {
            thisu1mEA = rowSums(sweep(probs, 2, STATS = 0:n, "*"))
            warning("Using median predictions is not supported for discrete distributions. Used mean prediction instead")
          }
          thisvar.estEA = rowSums(sweep(probs, 2, STATS = (0:n)^2, "*")) - thisu1mEA^2
          intervals = apply(probs, 1, generateBinomialInterval)
          thislowerEA = sapply(intervals, function(x) {x$lower})
          thisupperEA = sapply(intervals, function(x) {x$upper})
          thisRightRejectEA = sapply(intervals, function(x) {x$rightRejectProb})
          thisLeftRejectEA = sapply(intervals, function(x) {x$leftRejectProb})
        }
      }
      else {
        thisu1mEA = NA
        thislowerEA = NA
        thisupperEA = NA
        thisvar.estEA = NA
        thisRightRejectEA = NA
        thisLeftRejectEA = NA
      }
      
      # pixel level
      if(keepPixelPreds) {
        if(!genCountLevel) {
          if(predictionType == "mean")
            thisu1mPixel = rowMeans(pixelPredMat)
          else
            thisu1mPixel = apply(pixelPredMat, 1, median)
          thislowerPixel = logit(apply(pixelPredMat, 1, function(x) {quantile(x, probs=.1)}))
          thisupperPixel = logit(apply(pixelPredMat, 1, function(x) {quantile(x, probs=.9)}))
          thisvar.estPixel = apply(logit(pixelPredMat), 1, var)
          thisRightRejectPixel = NA
          thisLeftRejectPixel = NA
        } else {
          thisu1mPixel = pixelPredMat$preds
          thisvar.estPixel = pixelPredMat$vars
          thisupperPixel = pixelPredMat$upper
          thislowerPixel = pixelPredMat$lower
          thisRightRejectPixel = pixelPredMat$leftRejectProb
          thisLeftRejectPixel = pixelPredMat$rightRejectProb
        }
      }
      else {
        thisu1mPixel = NA
        thislowerPixel = NA
        thisupperPixel = NA
        thisvar.estPixel = NA
        thisRightRejectPixel = NA
        thisLeftRejectPixel = NA
      }
      
      # County level (county level results are always included)
      if(!genCountLevel) {
        if(predictionType == "mean")
          thisu1mCounty = rowMeans(countyPredMat)
        else
          thisu1mCounty = apply(countyPredMat, 1, median)
        thislowerCounty = logit(apply(countyPredMat, 1, function(x) {quantile(x, probs=.1)}))
        thisupperCounty = logit(apply(countyPredMat, 1, function(x) {quantile(x, probs=.9)}))
        thisvar.estCounty = apply(logit(countyPredMat), 1, var)
        thisRightRejectCounty = NA
        thisLeftRejectCounty = NA
      } else {
        thisu1mCounty = countyPredMat$preds
        thisvar.estCounty = countyPredMat$vars
        thisupperCounty = countyPredMat$upper
        thislowerCounty = countyPredMat$lower
        thisRightRejectCounty = countyPredMat$leftRejectProb
        thisLeftRejectCounty = countyPredMat$rightRejectProb
      }
      
      # region level
      if(genRegionLevel) {
        if(!genCountLevel) {
          if(predictionType == "mean")
            thisu1mRegion = rowMeans(regionPredMat)
          else
            thisu1mRegion = apply(regionPredMat, 1, median)
          thislowerRegion = logit(apply(regionPredMat, 1, function(x) {quantile(x, probs=.1)}))
          thisupperRegion = logit(apply(regionPredMat, 1, function(x) {quantile(x, probs=.9)}))
          thisvar.estRegion = apply(logit(regionPredMat), 1, var)
          thisRightRejectRegion = NA
          thisLeftRejectRegion = NA
        } else {
          thisu1mRegion = regionPredMat$preds
          thisvar.estRegion = regionPredMat$vars
          thisupperRegion = regionPredMat$upper
          thislowerRegion = regionPredMat$lower
          thisRightRejectRegion = regionPredMat$leftRejectProb
          thisLeftRejectRegion = regionPredMat$rightRejectProb
        }
      }
      else {
        thisu1mRegion = NA
        thislowerRegion = NA
        thisupperRegion = NA
        thisvar.estRegion = NA
        thisRightRejectRegion = NA
        thisLeftRejectRegion = NA
      }
      
      # collect results in a list, one element for each simulation
      countyResults[[i]] = data.frame(list(admin1=counties, 
                                           u1m.spde=thisu1mCounty, lower.spde=thislowerCounty, 
                                           upper.spde=thisupperCounty, logit.est.spde=logit(thisu1mCounty), 
                                           var.est.spde=thisvar.estCounty, 
                                           leftReject.spde=thisLeftRejectCounty, 
                                           rightReject.spde=thisRightRejectCounty))
      regionResults[[i]] = data.frame(list(region=regions, 
                                           u1m.spde=thisu1mRegion, lower.spde=thislowerRegion, 
                                           upper.spde=thisupperRegion, logit.est.spde=logit(thisu1mRegion), 
                                           var.est.spde=thisvar.estRegion, 
                                           leftReject.spde=thisLeftRejectRegion, 
                                           rightReject.spde=thisRightRejectRegion))
      pixelResults[[i]] = data.frame(list(u1m.spde=thisu1mPixel, lower.spde=thislowerPixel, 
                                          upper.spde=thisupperPixel, logit.est.spde=logit(thisu1mPixel), 
                                          var.est.spde=thisvar.estPixel, 
                                          leftReject.spde=thisLeftRejectPixel, 
                                          rightReject.spde=thisRightRejectPixel))
      eaResults[[i]] = data.frame(list(u1m.spde=thisu1mEA, lower.spde=thislowerEA, 
                                       upper.spde=thisupperEA, logit.est.spde=logit(thisu1mEA), 
                                       var.est.spde=thisvar.estEA, 
                                       leftReject.spde=thisLeftRejectEA, 
                                       rightReject.spde=thisRightRejectEA))
    }
  } else {
    ##### 
    ##### 
    ##### 
    ##### 
    ##### 
    # parallel version
    results = foreach(i = 1:nsim, .combine=rbind, .verbose=TRUE, .multicombine=TRUE) %dopar% {
      print(paste0("iteration ", i, "/", nsim))
      
      # get the simulated sample
      clustDat = clustDatMulti[[i]]
      
      # get observations from dataset
      obsCoords = cbind(clustDat$east, clustDat$north)
      obsNs = rep(25, nrow(obsCoords))
      obsCounts = clustDat$died
      obsUrban = clustDat$urban
      
      # fit model, get all predictions for each areal level and each posterior sample
      fit = fitSPDEModel(obsCoords, obsNs=obsNs, obsCounts, obsUrban, predCoords, predNs=predNs, 
                         predUrban, genCountyLevel=TRUE, popGrid=popGrid, nPostSamples=nPostSamples, 
                         verbose = verbose, includeClustEffect=includeClustEffect, 
                         int.strategy=int.strategy, genRegionLevel=genRegionLevel, 
                         keepPixelPreds=keepPixelPreds, genEALevel=genEALevel, 
                         urbanEffect=urbanEffect, link=1, predictionType=predictionType, 
                         exactAggregation=exactAggregation, genCountLevel=genCountLevel, 
                         eaDat=eaDat)
      print(paste0("Fit completed: iteration ", i, "/", nsim))
      countyPredMat = fit$countyPredMat
      regionPredMat = fit$regionPredMat
      pixelPredMat = fit$pixelPredMat
      eaPredMat = fit$eaPredMat
      eaMarginals = fit$eaMarginals
      
      ## calculate model fit properties for each level of predictions
      
      # enumeration area level
      if(genEALevel) {
        if(!genCountLevel) {
          if(predictionType == "mean")
            thisu1mEA = rowMeans(eaPredMat)
          else
            thisu1mEA = apply(eaPredMat, 1, median)
          thislowerEA = logit(apply(eaPredMat, 1, function(x) {quantile(x, probs=.1)}))
          thisupperEA = logit(apply(eaPredMat, 1, function(x) {quantile(x, probs=.9)}))
          thisvar.estEA = apply(logit(eaPredMat), 1, var)
          thisRightRejectEA = NA
          thisLeftRejectEA = NA
        }
        else {
          # get the marginal "binomial" densities at each location
          allMarginals = fit$mod$marginals.linear.predictor[fit$predInds]
          n = 25
          binomProb = function(k, i) {
            inla.emarginal(function(logitP) {dbinom(k, n, expit(logitP))}, allMarginals[[i]])
          }
          
          ## make highest density coverage interval on count scale
          # generate the "binomial" pmfs for each marginal
          probs = matrix(mapply(binomProb, 
                                k = matrix(0:n, nrow=length(allMarginals), ncol=n + 1, byrow = TRUE), 
                                i = matrix(1:length(allMarginals), nrow=length(allMarginals), ncol=n + 1)), 
                         nrow=length(allMarginals), ncol=n + 1)
          
          # now generate the summary statistics about the marginals
          if(predictionType == "mean")
            thisu1mEA = rowSums(sweep(probs, 2, STATS = 0:n, "*"))
          else {
            thisu1mEA = rowSums(sweep(probs, 2, STATS = 0:n, "*"))
            warning("Using median predictions is not supported for discrete distributions. Used mean prediction instead")
          }
          thisvar.estEA = rowSums(sweep(probs, 2, STATS = (0:n)^2, "*")) - thisu1mEA^2
          intervals = apply(probs, 1, generateBinomialInterval)
          thislowerEA = sapply(intervals, function(x) {x$lower})
          thisupperEA = sapply(intervals, function(x) {x$upper})
          thisRightRejectEA = sapply(intervals, function(x) {x$rightRejectProb})
          thisLeftRejectEA = sapply(intervals, function(x) {x$leftRejectProb})
        }
      }
      else {
        thisu1mEA = NA
        thislowerEA = NA
        thisupperEA = NA
        thisvar.estEA = NA
        thisRightRejectEA = NA
        thisLeftRejectEA = NA
      }
      print(paste0("EA results generated: iteration ", i, "/", nsim))
      
      # pixel level
      if(keepPixelPreds) {
        if(!genCountLevel) {
          if(predictionType == "mean")
            thisu1mPixel = rowMeans(pixelPredMat)
          else
            thisu1mPixel = apply(pixelPredMat, 1, median)
          thislowerPixel = logit(apply(pixelPredMat, 1, function(x) {quantile(x, probs=.1)}))
          thisupperPixel = logit(apply(pixelPredMat, 1, function(x) {quantile(x, probs=.9)}))
          thisvar.estPixel = apply(logit(pixelPredMat), 1, var)
          thisRightRejectPixel = NA
          thisLeftRejectPixel = NA
        } else {
          thisu1mPixel = pixelPredMat$preds
          thisvar.estPixel = pixelPredMat$vars
          thisupperPixel = pixelPredMat$upper
          thislowerPixel = pixelPredMat$lower
          thisRightRejectPixel = pixelPredMat$leftRejectProb
          thisLeftRejectPixel = pixelPredMat$rightRejectProb
        }
      }
      else {
        thisu1mPixel = NA
        thislowerPixel = NA
        thisupperPixel = NA
        thisvar.estPixel = NA
        thisRightRejectPixel = NA
        thisLeftRejectPixel = NA
      }
      print(paste0("Pixel results generated: iteration ", i, "/", nsim))
      
      # County level (county level results are always included)
      if(!genCountLevel) {
        if(predictionType == "mean")
          thisu1mCounty = rowMeans(countyPredMat)
        else
          thisu1mCounty = apply(countyPredMat, 1, median)
        thislowerCounty = logit(apply(countyPredMat, 1, function(x) {quantile(x, probs=.1)}))
        thisupperCounty = logit(apply(countyPredMat, 1, function(x) {quantile(x, probs=.9)}))
        thisvar.estCounty = apply(logit(countyPredMat), 1, var)
        thisRightRejectCounty = NA
        thisLeftRejectCounty = NA
      } else {
        thisu1mCounty = countyPredMat$preds
        thisvar.estCounty = countyPredMat$vars
        thisupperCounty = countyPredMat$upper
        thislowerCounty = countyPredMat$lower
        thisRightRejectCounty = countyPredMat$leftRejectProb
        thisLeftRejectCounty = countyPredMat$rightRejectProb
      }
      print(paste0("County results generated: iteration ", i, "/", nsim))
      
      # region level
      if(genRegionLevel) {
        if(!genCountLevel) {
          if(predictionType == "mean")
            thisu1mRegion = rowMeans(regionPredMat)
          else
            thisu1mRegion = apply(regionPredMat, 1, median)
          thislowerRegion = logit(apply(regionPredMat, 1, function(x) {quantile(x, probs=.1)}))
          thisupperRegion = logit(apply(regionPredMat, 1, function(x) {quantile(x, probs=.9)}))
          thisvar.estRegion = apply(logit(regionPredMat), 1, var)
          thisRightRejectRegion = NA
          thisLeftRejectRegion = NA
        } else {
          thisu1mRegion = regionPredMat$preds
          thisvar.estRegion = regionPredMat$vars
          thisupperRegion = regionPredMat$upper
          thislowerRegion = regionPredMat$lower
          thisRightRejectRegion = regionPredMat$leftRejectProb
          thisLeftRejectRegion = regionPredMat$rightRejectProb
        }
      }
      else {
        thisu1mRegion = NA
        thislowerRegion = NA
        thisupperRegion = NA
        thisvar.estRegion = NA
        thisRightRejectRegion = NA
        thisLeftRejectRegion = NA
      }
      print(paste0("Region results generated: iteration ", i, "/", nsim))
      
      # collect results in a list, one element for each simulation
      thisCountyResults = data.frame(list(admin1=counties, 
                                          u1m.spde=thisu1mCounty, lower.spde=thislowerCounty, 
                                          upper.spde=thisupperCounty, logit.est.spde=logit(thisu1mCounty), 
                                          var.est.spde=thisvar.estCounty, 
                                          leftReject.spde=thisLeftRejectCounty, 
                                          rightReject.spde=thisRightRejectCounty))
      thisRegionResults = data.frame(list(region=regions, 
                                          u1m.spde=thisu1mRegion, lower.spde=thislowerRegion, 
                                          upper.spde=thisupperRegion, logit.est.spde=logit(thisu1mRegion), 
                                          var.est.spde=thisvar.estRegion, 
                                          leftReject.spde=thisLeftRejectRegion, 
                                          rightReject.spde=thisRightRejectRegion))
      thisPixelResults = data.frame(list(u1m.spde=thisu1mPixel, lower.spde=thislowerPixel, 
                                         upper.spde=thisupperPixel, logit.est.spde=logit(thisu1mPixel), 
                                         var.est.spde=thisvar.estPixel, 
                                         leftReject.spde=thisLeftRejectPixel, 
                                         rightReject.spde=thisRightRejectPixel))
      thisEaResults = data.frame(list(u1m.spde=thisu1mEA, lower.spde=thislowerEA, 
                                      upper.spde=thisupperEA, logit.est.spde=logit(thisu1mEA), 
                                      var.est.spde=thisvar.estEA, 
                                      leftReject.spde=thisLeftRejectEA, 
                                      rightReject.spde=thisRightRejectEA))
      cbind(thisCountyResults, thisRegionResults, thisPixelResults, thisEaResults)
    }
    
    # separate results into the different aggregation levels
    countyResults = results[,1:8]
    regionResults = results[,9:16]
    pixelResults = results[,17:23]
    eaResults = results[,24:30]
  }
  list(countyResults=countyResults, regionResults=regionResults, pixelResults=pixelResults, 
       eaResults=eaResults)
}

# compute true proportion of population that died in each county in the order of the 
# input counties variable from the full dataset eaDat
getTruthByCounty = function(eaDat, counties=as.character(unique(mort$admin1))) {
  # get true mortality rate for each county
  # diedByCounty = aggregate(eaDat$died, by=list(eaDat$admin1), FUN="sum")
  # theseCounties = diedByCounty$Group.1
  # diedByCounty = diedByCounty$V1
  diedByCounty = tapply(eaDat$died, eaDat$admin1, sum)
  theseCounties = rownames(diedByCounty)
  # nByCounty = aggregate(eaDat$numChildren, by=list(eaDat$admin1), FUN="sum")$x
  nByCounty = tapply(eaDat$numChildren, eaDat$admin1, sum)
  mortRate = diedByCounty/nByCounty
  
  # sort to be in the order of the input countries variable
  sortI = match(counties, theseCounties)
  list(counties=counties, mortRate=mortRate[sortI])
}