# code for generating table for Jon's talk with model scores

# generate predictive results given the simulated datasets simDataMulti.RData
# includeClustEffect, urbanEffect: control whether or knot to include a 
# cluster random effect in urban fixed effect in the model
# genRegionLevel, keepPixelPreds, genEALevel: control whether to include 
# regional, pixel, and enumeration area level predictions
resultsSPDE = function(nPostSamples=100, test=FALSE, nTest=5, verbose=TRUE, 
                       includeClustEffect=FALSE, int.strategy="eb", 
                       genRegionLevel=TRUE, keepPixelPreds=TRUE, 
                       genEALevel=TRUE, urbanEffect=TRUE) {
  # Load data
  # load("simDataMulti.RData") # overSampDat, SRSDat
  # load a different 1 of these depending on whether a cluster effect should be included 
  # in the simulation of the data or not (tausq is the cluster effect variance)
  load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOver2.RData")
  # load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOver2.RData")
  eaDat = overSampDat$eaDat
  clustSRS = SRSDat$clustDat
  clustOverSamp = overSampDat$clustDat
  
  if(test) {
    clustSRS = lapply(1:nTest, function(i) {clustSRS[[i]]})
    clustOverSamp = lapply(1:nTest, function(i) {clustOverSamp[[i]]})
  }
  
  spdeSRS = resultsSPDEHelper(clustSRS, eaDat, nPostSamples = nPostSamples, verbose=verbose, 
                              includeClustEffect=includeClustEffect, int.strategy=int.strategy, 
                              genRegionLevel=genRegionLevel, keepPixelPreds=keepPixelPreds, 
                              genEALevel=genEALevel, urbanEffect=urbanEffect)
  spdeOverSamp = resultsSPDEHelper(clustSRS, eaDat, nPostSamples = nPostSamples, verbose=verbose, 
                                   includeClustEffect=includeClustEffect, int.strategy=int.strategy, 
                                   genRegionLevel=genRegionLevel, keepPixelPreds=keepPixelPreds, 
                                   genEALevel=genEALevel, urbanEffect=urbanEffect)
  
  # save(spdeSRS, spdeOverSamp, file="resultsSPDETausq0.RData")
  # save(spdeSRS, spdeOverSamp, file="resultsSPDETausq0.01.RData")
  
  list(spdeSRS=spdeSRS, spdeOverSamp=spdeOverSamp)
}

# generate predictive results given a simulated dataset simDataMulti.RData
resultsSPDEHelper = function(clustDatMulti, eaDat, nPostSamples=100, verbose=TRUE, 
                             includeClustEffect=FALSE, int.strategy="eb", 
                             genRegionLevel=TRUE, keepPixelPreds=TRUE, genEALevel=TRUE, 
                             urbanEffect=TRUE) {
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
  popGrid = makeInterpPopGrid()
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
  countyResults = list()
  regionResults = list()
  pixelResults = list()
  eaResults = list()
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
                       urbanEffect=urbanEffect)
    countyPredMat = fit$countyPredMat
    regionPredMat = fit$regionPredMat
    pixelPredMat = fit$pixelPredMat
    eaPredMat = fit$eaPredMat
    
    ## calculate model fit properties for each level of predictions
    # County level (this is always included)
    thisu1mCounty = rowMeans(countyPredMat)
    thislowerCounty = logit(apply(countyPredMat, 1, function(x) {quantile(x, probs=.1)}))
    thisupperCounty = logit(apply(countyPredMat, 1, function(x) {quantile(x, probs=.9)}))
    thisvar.estCounty = apply(logit(countyPredMat), 1, var)
    
    # region level
    if(genRegionLevel) {
      thisu1mRegion = rowMeans(regionPredMat)
      thislowerRegion = logit(apply(regionPredMat, 1, function(x) {quantile(x, probs=.1)}))
      thisupperRegion = logit(apply(regionPredMat, 1, function(x) {quantile(x, probs=.9)}))
      thisvar.estRegion = apply(logit(regionPredMat), 1, var)
    }
    else {
      thisu1mRegion = NA
      thislowerRegion = NA
      thisupperRegion = NA
      thisvar.estRegion = NA
    }
    
    # pixel level
    if(keepPixelPreds) {
      thisu1mPixel = rowMeans(pixelPredMat)
      thislowerPixel = logit(apply(pixelPredMat, 1, function(x) {quantile(x, probs=.1)}))
      thisupperPixel = logit(apply(pixelPredMat, 1, function(x) {quantile(x, probs=.9)}))
      thisvar.estPixel = apply(logit(pixelPredMat), 1, var)
    }
    else {
      thisu1mPixel = NA
      thislowerPixel = NA
      thisupperPixel = NA
      thisvar.estPixel = NA
    }
    
    # enumeration area level
    if(genEALevel) {
      thisu1mEA = rowMeans(eaPredMat)
      thislowerEA = logit(apply(eaPredMat, 1, function(x) {quantile(x, probs=.1)}))
      thisupperEA = logit(apply(eaPredMat, 1, function(x) {quantile(x, probs=.9)}))
      thisvar.estEA = apply(logit(eaPredMat), 1, var)
    }
    else {
      thisu1mEA = NA
      thislowerEA = NA
      thisupperEA = NA
      thisvar.estEA = NA
    }
    
    # collect results in a list, one element for each simulation
    countyResults = c(countyResults, 
                      list(data.frame(list(admin1=counties, 
                                           u1m.spde=thisu1mCounty, lower.spde=thislowerCounty, 
                                           upper.spde=thisupperCounty, logit.est.spde=logit(thisu1mCounty), 
                                           var.est.spde=thisvar.estCounty))))
    regionResults = c(regionResults, 
                      list(data.frame(list(region=regions, 
                                           u1m.spde=thisu1mRegion, lower.spde=thislowerRegion, 
                                           upper.spde=thisupperRegion, logit.est.spde=logit(thisu1mRegion), 
                                           var.est.spde=thisvar.estRegion))))
    pixelResults = c(pixelResults, 
                      list(data.frame(list(u1m.spde=thisu1mPixel, lower.spde=thislowerPixel, 
                                           upper.spde=thisupperPixel, logit.est.spde=logit(thisu1mPixel), 
                                           var.est.spde=thisvar.estPixel))))
    eaResults = c(eaResults, 
                      list(data.frame(list(u1m.spde=thisu1mEA, lower.spde=thislowerEA, 
                                           upper.spde=thisupperEA, logit.est.spde=logit(thisu1mEA), 
                                           var.est.spde=thisvar.estEA))))
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