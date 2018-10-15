##### pick the following settings before running code:
## pick cluster variance and load associated superpopulation:
# EITHER: cluster variance = 0.01
tausq = .1^2
out = load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOver2.RData")
eaDat = SRSDat$eaDat

# OR: cluster variance = 0 (no cluster effect)
tausq = 0
out = load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOver2.RData")
eaDat = SRSDat$eaDat

## pick one of the subLevels to get results at:
resultType = "county"
resultType = "pixel" # not yet working
resultType = "EA"

## pick which sampling mechanism:
# EITHER: SRS
sampling = "SRS"
clustDat = SRSDat

# OR: urban oversampling
sampling = "oversamp"
clustDat = overSampDat

## generate the true values at the county level for the 
## given settings if these are new settings
if(0){
  if(resultType == "county") {
    # compute truth based on superpopulation
    regions = sort(unique(eaDat$admin1))
    truthbycounty <- rep(NA, 47)
    
    for(i in 1:47){
      super = eaDat[eaDat$admin1 == regions[i],]
      truthbycounty[i] <- sum(super$died)/sum(super$numChildren)
    }
    truth = data.frame(admin1=regions, truth=truthbycounty)
    save(truth, file="truthbycounty.RData")
  } else if(resultType == "pixel") {
    counties = sort(unique(eaDat$admin1))
    eaToPixel = eaDat$pixelI
    childrenPerPixel = tapply(eaDat$numChildren, list(pixel=eaDat$pixelI), sum)
    urbanPixel = tapply(eaDat$urban, list(pixel=eaDat$pixelI), function(x) {mean(x[1])})
    deathsPerPixel = tapply(eaDat$died, list(pixel=eaDat$pixelI), sum)
    regions = names(childrenPerPixel) # these are the pixels with enumeration areas in them
    pixelToAdmin = match(popGrid$admin1[as.numeric(regions)], counties)
    
    truth = data.frame(pixel=regions, truth=deathsPerPixel / childrenPerPixel, countyI=pixelToAdmin, urban=urbanPixel)
    save(truth, file="truthbyPixel.RData")
  } else if(resultType == "EA") {
    truth = data.frame(EA = 1:nrow(eaDat), truth = eaDat$died/25, urban=eaDat$urban)
    save(truth, file="truthbyEA.RData")
  }
}

##### Now run the rest of the script (except possibly for the plotting, depending 
##### on the result aggregation level).

# load data
if(tausq == .1^2) {
  out = load("resultsDirectNaiveTausq0.01.RData")
  out = load("resultsMercerTausq0.01.RData")
  out = load("KenyaSpatialDesignResultNewTausq0.01UrbRurFALSE.RData")
  designResNoUrb = designRes
  out = load("KenyaSpatialDesignResultNewTausq0.01UrbRurTRUE.RData")
  out = load("resultsSPDETausq0.01urbanEffectFALSE.RData")
  spdeSRSNoUrb = spdeSRS
  spdeOverSampNoUrb = spdeOverSamp
  out = load("resultsSPDETausq0.01urbanEffectTRUE.RData")
} else if(tausq == 0) {
  out = load("resultsDirectNaiveTausq0.RData")
  out = load("resultsMercerTausq0.RData")
  out = load("KenyaSpatialDesignResultNewTausq0UrbRurFALSE.RData")
  designResNoUrb = designRes
  out = load("KenyaSpatialDesignResultNewTausq0UrbRurTRUE.RData")
  out = load("resultsSPDETausq0urbanEffectFALSE.RData")
  spdeSRSNoUrb = spdeSRS
  spdeOverSampNoUrb = spdeOverSamp
  out = load("resultsSPDETausq0urbanEffectTRUE.RData")
}

# modify discrete estimates stratified by urban/rural if resultType is at finer scale than county
# (currently, further resolving discrete estimates is only supported for the direct estimates)
if(resultType == "EA") {
  filterByUrban = function(i) {
    urbanI = 1:nrow(directEstSRS[[i]])
  }
  directEstSRS = lapply(1:100, )
}

source("scores.R")

# load the truth
load(paste0("truthby", resultType, ".RData"))

# function for generating discrete model pixel and EA level predictions 
# given the county level predictions
getSubLevelResults = function(resultTable) {
  # get the order of the counties that resultTable is in
  counties = sort(unique(eaDat$admin1))
  
  # convert results to the desired aggregation level if necessary for discrete models:
  if(resultType == "pixel") {
    # compute pixel results
    pixelToAdmin = match(popGrid$admin1, counties)
    resultTable = resultTable[pixelToAdmin,]
  }
  else if(resultType == "EA") {
    # compute EA level results
    eaToAdmin = match(eaDat$admin1, counties)
    resultTable = resultTable[eaToAdmin,]
  }
  
  # modify result row and column table names according to aggregation level
  whichName = which(names(resultTable) == "admin1")
  names(resultTable)[whichName] = resultType
  if((resultType == "EA") && (is.data.frame(resultTable))) {
    # resultTable[[resultType]] = resultTable$eaI
    resultTable[[resultType]] = 1:nrow(resultTable)
  } else if((resultType == "pixel") && (is.data.frame(resultTable))) {
    # resultTable[[resultType]] = resultTable$eaI
    resultTable[[resultType]] = 1:nrow(resultTable)
  }
  
  if(resultType == "pixel")
    resultTable = resultTable[as.numeric(as.character(truth$pixel)),]
  
  resultTable
}

# convert the truth to the desired aggregation level
if(resultType == "county")
  truth = getSubLevelResults(truth)

# compute scores
scoresDirectSRS = scoresNaiveSRS = scoresMercerSRS = scoresBYMNoUrbSRS = scoresBYMSRS = scoresSPDENoUrbSRS = scoresSPDESRS = data.frame()
scoresDirectoverSamp = scoresNaiveoverSamp = scoresMerceroverSamp = scoresBYMNoUrboverSamp = scoresBYMoverSamp = scoresSPDENoUrboverSamp = scoresSPDEoverSamp = data.frame()
# scoresDirectSRS = scoresNaiveSRS = scoresMercerSRS = scoresBYMSRS = data.frame()
# scoresDirectoverSamp = scoresNaiveoverSamp = scoresMerceroverSamp = scoresBYMoverSamp = data.frame()

# convert results to the desired aggregation level
# not including urban effect
designResNoUrb$SRSdat$Q10 = getSubLevelResults(designResNoUrb$SRSdat$Q10)
designResNoUrb$SRSdat$Q50 = getSubLevelResults(designResNoUrb$SRSdat$Q50)
designResNoUrb$SRSdat$Q90 = getSubLevelResults(designResNoUrb$SRSdat$Q90)
designResNoUrb$SRSdat$mean = getSubLevelResults(designResNoUrb$SRSdat$mean)
designResNoUrb$SRSdat$stddev = getSubLevelResults(designResNoUrb$SRSdat$stddev)
designResNoUrb$overSampDat$Q10 = getSubLevelResults(designResNoUrb$overSampDat$Q10)
designResNoUrb$overSampDat$Q50 = getSubLevelResults(designResNoUrb$overSampDat$Q50)
designResNoUrb$overSampDat$Q90 = getSubLevelResults(designResNoUrb$overSampDat$Q90)
designResNoUrb$overSampDat$mean = getSubLevelResults(designResNoUrb$overSampDat$mean)
designResNoUrb$overSampDat$stddev = getSubLevelResults(designResNoUrb$overSampDat$stddev)

# including urban effect
designRes$SRSdat$Q10 = getSubLevelResults(designRes$SRSdat$Q10)
designRes$SRSdat$Q50 = getSubLevelResults(designRes$SRSdat$Q50)
designRes$SRSdat$Q90 = getSubLevelResults(designRes$SRSdat$Q90)
designRes$SRSdat$mean = getSubLevelResults(designRes$SRSdat$mean)
designRes$SRSdat$stddev = getSubLevelResults(designRes$SRSdat$stddev)
designRes$overSampDat$Q10 = getSubLevelResults(designRes$overSampDat$Q10)
designRes$overSampDat$Q50 = getSubLevelResults(designRes$overSampDat$Q50)
designRes$overSampDat$Q90 = getSubLevelResults(designRes$overSampDat$Q90)
designRes$overSampDat$mean = getSubLevelResults(designRes$overSampDat$mean)
designRes$overSampDat$stddev = getSubLevelResults(designRes$overSampDat$stddev)

for(i in c(1:100)) { # for problem fitting mercerSRS for SRS sampling, tausq=0
  # for(i in 1:100) {
  print(i)
  resultName = paste0(resultType, "Results")
  if(resultType == "EA")
    resultName = "eaResults"
  
  # convert results to the desired aggregation level
  if(sampling == "SRS") {
    directEstSRSi = getSubLevelResults(directEstSRS[[i]])
    naiveSRSi = getSubLevelResults(naiveSRS[[i]])
    mercerSRSi = getSubLevelResults(mercerSRS[[i]])
    if(resultType != "county") {
      spdeSRSNoUrbi = spdeSRSNoUrb[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
      spdeSRSi = spdeSRS[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
      directEstSRSUrbani = getSubLevelResults(directEstSRSUrb[[i]])
      directEstSRSRurali = getSubLevelResults(directEstSRSRural[[i]])
    } else {
      spdeSRSNoUrbi = spdeSRSNoUrb[[resultName]][[i]]
      spdeSRSi = spdeSRS[[resultName]][[i]]
    }
  } else if(sampling == "oversamp") {
    directEstoverSampi = getSubLevelResults(directEstoverSamp[[i]])
    naiveoverSampi = getSubLevelResults(naiveoverSamp[[i]])
    merceroverSampi = getSubLevelResults(merceroverSamp[[i]])
    if(resultType != "county") {
      spdeOverSampNoUrbi = spdeOverSampNoUrb[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
      spdeOverSampi = spdeOverSamp[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
      directEstoverSampUrbani = getSubLevelResults(directEstoverSampUrb[[i]])
      directEstoverSampRurali = getSubLevelResults(directEstoverSampRural[[i]])
    } else {
      spdeOverSampNoUrbi = spdeOverSampNoUrb[[resultName]][[i]]
      spdeOverSampi = spdeOverSamp[[resultName]][[i]]
    }
  }
  
  if(resultType == "EA") {
    # set first row of spde results to be the EA index
    if(sampling == "SRS") {
      spdeSRSNoUrbi[[resultType]] = 1:nrow(spdeSRSNoUrbi)
      spdeSRSi[[resultType]] = 1:nrow(spdeSRSi)
      
      whichName = which(names(spdeSRSNoUrbi) == "EA")
      spdeSRSNoUrbi = cbind(spdeSRSNoUrbi[,whichName], spdeSRSNoUrbi[,-whichName])
      whichName = which(names(spdeSRSi) == "EA")
      spdeSRSi = cbind(spdeSRSi[,whichName], spdeSRSi[,-whichName])
      
      names(spdeSRSNoUrbi)[1] = "EA"
      names(spdeSRSi)[1] = "EA"
    } else {
      spdeOverSampNoUrbi[[resultType]] = 1:nrow(spdeOverSampNoUrbi)
      spdeOverSampi[[resultType]] = 1:nrow(spdeOverSampi)
      
      whichName = which(names(spdeOverSampNoUrbi) == "EA")
      spdeOverSampNoUrbi = cbind(spdeOverSampNoUrbi[,whichName], spdeOverSampNoUrbi[,-whichName])
      whichName = which(names(spdeOverSampi) == "EA")
      spdeOverSampi = cbind(spdeOverSampi[,whichName], spdeOverSampi[,-whichName])
      
      names(spdeOverSampNoUrbi)[1] = "EA"
      names(spdeOverSampi)[1] = "EA"
    }
  }
  
  # for spde results, modify the name of the results
  # modify result row and column table names according to aggregation level
  if(resultType == "county") {
    # without urban effect:
    if(sampling == "SRS") {
      whichName = which(names(spdeSRSNoUrbi) == "admin1")
      names(spdeSRSNoUrbi)[whichName] = resultType
    } else if (sampling == "oversamp") {
      whichName = which(names(spdeOverSampNoUrbi) == "admin1")
      names(spdeOverSampNoUrbi)[whichName] = resultType
    }
    
    # with urban effect:
    if(sampling == "SRS") {
      whichName = which(names(spdeSRSi) == "admin1")
      names(spdeSRSi)[whichName] = resultType
    } else if(sampling == "oversamp") {
      whichName = which(names(spdeOverSampi) == "admin1")
      names(spdeOverSampi)[whichName] = resultType
    }
  }
  
  if(resultType == "pixel") {
    # set first row of spde results to be the pixel index
    if(sampling == "SRS") {
      spdeSRSNoUrbi[[resultType]] = truth$pixel
      spdeSRSi[[resultType]] = truth$pixel
      
      whichName = which(names(spdeSRSNoUrbi) == "pixel")
      spdeSRSNoUrbi = cbind(spdeSRSNoUrbi[,whichName], spdeSRSNoUrbi[,-whichName])
      whichName = which(names(spdeSRSi) == "pixel")
      spdeSRSi = cbind(spdeSRSi[,whichName], spdeSRSi[,-whichName])
      
      names(spdeSRSNoUrbi)[1] = "pixel"
      names(spdeSRSi)[1] = "pixel"
    } else {
      spdeOverSampNoUrbi[[resultType]] = truth$pixel
      spdeOverSampi[[resultType]] = truth$pixel
      
      whichName = which(names(spdeOverSampNoUrbi) == "pixel")
      spdeOverSampNoUrbi = cbind(spdeOverSampNoUrbi[,whichName], spdeOverSampNoUrbi[,-whichName])
      whichName = which(names(spdeOverSampi) == "pixel")
      spdeOverSampi = cbind(spdeOverSampi[,whichName], spdeOverSampi[,-whichName])
      
      names(spdeOverSampNoUrbi)[1] = "pixel"
      names(spdeOverSampi)[1] = "pixel"
    }
  }
  
  # change names of table variables in spde model with no urban effect to reflect that
  if(sampling == "SRS") {
    names(spdeSRSNoUrbi)[2:6] = paste0(names(spdeSRSNoUrbi)[2:6], "NoUrb")
  } else {
    names(spdeOverSampNoUrbi)[2:6] = paste0(names(spdeOverSampNoUrbi)[2:6], "NoUrb")
  }
  
  if(sampling == "SRS") {
    allresSRS = merge(truth, directEstSRSi, by=resultType)
    colnames(allresSRS) = c(resultType, "truth", paste(colnames(allresSRS)[3:8], "direct", sep=""))
    allresSRS = merge(allresSRS, naiveSRSi, by=resultType)
    allresSRS = merge(allresSRS, mercerSRSi, by=resultType)
    allresSRS = merge(allresSRS, spdeSRSNoUrbi, by=resultType)
    allresSRS = merge(allresSRS, spdeSRSi, by=resultType)
  }
  
  if(sampling == "oversamp") {
    allresoverSamp = merge(truth, directEstoverSampi, by=resultType)
    colnames(allresoverSamp) = c(resultType, "truth", paste(colnames(allresoverSamp)[3:8], "direct", sep=""))
    allresoverSamp = merge(allresoverSamp, naiveoverSampi, by=resultType)
    allresoverSamp = merge(allresoverSamp, merceroverSampi, by=resultType)
    allresoverSamp = merge(allresoverSamp, spdeOverSampNoUrbi, by=resultType)
    allresoverSamp = merge(allresoverSamp, spdeOverSampi, by=resultType)
  }
  
  # set whether or not to calculate scores on logit scale depending on result type
  thisTruth = allresSRS$truth
  useLogit = FALSE
  if(resultType != "EA" && resultType != "pixel") {
    thisTruth = logit(thisTruth)
    useLogit=TRUE
  }
  
  # SRS setting 
  if(sampling == "SRS") {
    my.biasSRSdirect = bias(thisTruth, allresSRS$logit.estdirect, logit=useLogit, my.var=allresSRS$var.estdirect)
    my.mseSRSdirect = mse(thisTruth, allresSRS$logit.estdirect, logit=useLogit, my.var=allresSRS$var.estdirect)
    my.dssSRSdirect = dss(thisTruth, allresSRS$logit.estdirect, allresSRS$var.estdirect)
    my.crpsSRSdirect = crpsNormal(thisTruth, allresSRS$logit.estdirect, allresSRS$var.estdirect, resultType=resultType)
    my.coverageSRSdirect = coverage(thisTruth, allresSRS$lowerdirect, allresSRS$upperdirect, logit=useLogit)
    my.lengthSRSdirect = intervalWidth(allresSRS$lowerdirect, allresSRS$upperdirect, logit=useLogit)
    
    my.biasSRSnaive = bias(thisTruth, allresSRS$logit.est, logit=useLogit, my.var=allresSRS$var.est)
    my.mseSRSnaive = mse(thisTruth, allresSRS$logit.est, logit=useLogit, my.var=allresSRS$var.est)
    my.dssSRSnaive = dss(thisTruth, allresSRS$logit.est, allresSRS$var.est)
    my.crpsSRSnaive = crpsNormal(thisTruth, allresSRS$logit.est, allresSRS$var.est, resultType=resultType)
    my.coverageSRSnaive = coverage(thisTruth, allresSRS$lower, allresSRS$upper, logit=useLogit)
    my.lengthSRSnaive = intervalWidth(allresSRS$lower, allresSRS$upper, logit=useLogit)
    
    my.biasSRSmercer = bias(thisTruth, allresSRS$logit.est.mercer, logit=useLogit, my.var=allresSRS$var.est.mercer)
    my.mseSRSmercer = mse(thisTruth, allresSRS$logit.est.mercer, logit=useLogit, my.var=allresSRS$var.est.mercer)
    my.dssSRSmercer = dss(thisTruth, allresSRS$logit.est.mercer, allresSRS$var.est.mercer)
    my.crpsSRSmercer = crpsNormal(thisTruth, allresSRS$logit.est.mercer, allresSRS$var.est.mercer, resultType=resultType)
    my.coverageSRSmercer = coverage(thisTruth, allresSRS$lower.mercer, allresSRS$upper.mercer, logit=useLogit)
    my.lengthSRSmercer = intervalWidth(allresSRS$lower.mercer, allresSRS$upper.mercer, logit=useLogit)
    
    my.biasSRSbymNoUrb = bias(thisTruth, designResNoUrb$SRSdat$mean[,i], logit=useLogit, my.var=(designResNoUrb$SRSdat$stddev[,i])^2)
    my.mseSRSbymNoUrb = mse(thisTruth, designResNoUrb$SRSdat$mean[,i], logit=useLogit, my.var=(designResNoUrb$SRSdat$stddev[,i])^2)
    my.dssSRSbymNoUrb = dss(thisTruth, designResNoUrb$SRSdat$mean[,i], (designResNoUrb$SRSdat$stddev[,i])^2)
    my.crpsSRSbymNoUrb = crpsNormal(thisTruth, designResNoUrb$SRSdat$mean[,i], (designResNoUrb$SRSdat$stddev[,i])^2, resultType=resultType)
    my.coverageSRSbymNoUrb = coverage(thisTruth, designResNoUrb$SRSdat$Q10[,i],designResNoUrb$SRSdat$Q90[,i], logit=useLogit)
    my.lengthSRSbymNoUrb = intervalWidth(designResNoUrb$SRSdat$Q10[,i], designResNoUrb$SRSdat$Q90[,i], logit=useLogit)
    
    my.biasSRSbym = bias(thisTruth, designRes$SRSdat$mean[,i], logit=useLogit, my.var=(designRes$SRSdat$stddev[,i])^2)
    my.mseSRSbym = mse(thisTruth, designRes$SRSdat$mean[,i], logit=useLogit, my.var=(designRes$SRSdat$stddev[,i])^2)
    my.dssSRSbym = dss(thisTruth, designRes$SRSdat$mean[,i], (designRes$SRSdat$stddev[,i])^2)
    my.crpsSRSbym = crpsNormal(thisTruth, designRes$SRSdat$mean[,i], (designRes$SRSdat$stddev[,i])^2, resultType=resultType)
    my.coverageSRSbym = coverage(thisTruth, designRes$SRSdat$Q10[,i],designRes$SRSdat$Q90[,i], logit=useLogit)
    my.lengthSRSbym = intervalWidth(designRes$SRSdat$Q10[,i], designRes$SRSdat$Q90[,i], logit=useLogit)
    
    my.biasSRSspdeNoUrb = bias(thisTruth, allresSRS$logit.est.spdeNoUrb, logit=useLogit, my.var=allresSRS$var.est.spdeNoUrb)
    my.mseSRSspdeNoUrb = mse(thisTruth, allresSRS$logit.est.spdeNoUrb, logit=useLogit, my.var=allresSRS$var.est.spdeNoUrb)
    my.dssSRSspdeNoUrb = dss(thisTruth, allresSRS$logit.est.spdeNoUrb, allresSRS$var.est.spdeNoUrb)
    my.crpsSRSspdeNoUrb = crpsNormal(thisTruth, allresSRS$logit.est.spdeNoUrb, allresSRS$var.est.spdeNoUrb, resultType=resultType)
    my.coverageSRSspdeNoUrb = coverage(thisTruth, allresSRS$lower.spdeNoUrb, allresSRS$upper.spdeNoUrb, logit=useLogit)
    my.lengthSRSspdeNoUrb = intervalWidth(allresSRS$lower.spdeNoUrb, allresSRS$upper.spdeNoUrb, logit=useLogit)
    
    my.biasSRSspde = bias(thisTruth, allresSRS$logit.est.spde, logit=useLogit, my.var=allresSRS$var.est.spde)
    my.mseSRSspde = mse(thisTruth, allresSRS$logit.est.spde, logit=useLogit, my.var=allresSRS$var.est.spde)
    my.dssSRSspde = dss(thisTruth, allresSRS$logit.est.spde, allresSRS$var.est.spde)
    my.crpsSRSspde = crpsNormal(thisTruth, allresSRS$logit.est.spde, allresSRS$var.est.spde, resultType=resultType)
    my.coverageSRSspde = coverage(thisTruth, allresSRS$lower.spde, allresSRS$upper.spde, logit=useLogit)
    my.lengthSRSspde = intervalWidth(allresSRS$lower.spde, allresSRS$upper.spde, logit=useLogit)
    
    scoresDirectSRS <- rbind(scoresDirectSRS,
                             data.frame(dataset=i, region=allresSRS[[resultType]],
                                        bias=my.biasSRSdirect,
                                        mse=my.mseSRSdirect,
                                        dss=my.dssSRSdirect,
                                        coverage=my.coverageSRSdirect,
                                        var=mean(allresSRS$var.estdirect),
                                        crps=my.crpsSRSdirect, 
                                        length=my.lengthSRSdirect))
    
    scoresNaiveSRS <- rbind(scoresNaiveSRS,
                            data.frame(dataset=i, region=allresSRS[[resultType]],
                                       bias=my.biasSRSnaive,
                                       mse=my.mseSRSnaive,
                                       dss=my.dssSRSnaive,
                                       coverage=my.coverageSRSnaive,
                                       var=mean(allresSRS$var.est),
                                       crps=my.crpsSRSnaive, 
                                       length=my.lengthSRSnaive))
    scoresMercerSRS <- rbind(scoresMercerSRS,
                             data.frame(dataset=i, region=allresSRS[[resultType]],
                                        bias=my.biasSRSmercer,
                                        mse=my.mseSRSmercer,
                                        dss=my.dssSRSmercer,
                                        coverage=my.coverageSRSmercer,
                                        var=mean(allresSRS$var.est.mercer),
                                        crps=my.crpsSRSmercer, 
                                        length=my.lengthSRSmercer))
    scoresBYMNoUrbSRS <- rbind(scoresBYMNoUrbSRS,
                               data.frame(dataset=i, region=allresSRS[[resultType]],
                                          bias=my.biasSRSbymNoUrb,
                                          mse=my.mseSRSbymNoUrb,
                                          dss=my.dssSRSbymNoUrb,
                                          coverage=my.coverageSRSbymNoUrb,
                                          var=mean((designResNoUrb$SRSdat$stddev[,i])^2),
                                          crps=my.crpsSRSbymNoUrb, 
                                          length=my.lengthSRSbymNoUrb))
    scoresBYMSRS <- rbind(scoresBYMSRS,
                          data.frame(dataset=i, region=allresSRS[[resultType]],
                                     bias=my.biasSRSbym,
                                     mse=my.mseSRSbym,
                                     dss=my.dssSRSbym,
                                     coverage=my.coverageSRSbym,
                                     var=mean((designRes$SRSdat$stddev[,i])^2),
                                     crps=my.crpsSRSbym, 
                                     length=my.lengthSRSbym))
    scoresSPDENoUrbSRS <- rbind(scoresSPDENoUrbSRS, 
                                data.frame(dataset=i, region=allresSRS[[resultType]], 
                                           bias=my.biasSRSspdeNoUrb, 
                                           mse=my.mseSRSspdeNoUrb,
                                           dss=my.dssSRSspdeNoUrb,
                                           coverage=my.coverageSRSspdeNoUrb, 
                                           var=mean(allresSRS$var.est.spdeNoUrb),
                                           crps=my.crpsSRSspdeNoUrb, 
                                           length=my.lengthSRSspdeNoUrb))
    scoresSPDESRS <- rbind(scoresSPDESRS, 
                           data.frame(dataset=i, region=allresSRS[[resultType]], 
                                      bias=my.biasSRSspde, 
                                      mse=my.mseSRSspde,
                                      dss=my.dssSRSspde,
                                      coverage=my.coverageSRSspde, 
                                      var=mean(allresSRS$var.est.spde),
                                      crps=my.crpsSRSspde, 
                                      length=my.lengthSRSspde))
  } else if(sampling == "oversamp") {
    # oversampling setting
    my.biasoverSampdirect = bias(thisTruth, allresoverSamp$logit.estdirect, logit=useLogit, my.var=allresoverSamp$var.estdirect)
    my.mseoverSampdirect = mse(thisTruth, allresoverSamp$logit.estdirect, logit=useLogit, my.var=allresoverSamp$var.estdirect)
    my.dssoverSampdirect = dss(thisTruth, allresoverSamp$logit.estdirect, allresoverSamp$var.estdirect)
    my.crpsoverSampdirect = crpsNormal(thisTruth, allresoverSamp$logit.estdirect, allresoverSamp$var.estdirect, resultType=resultType)
    my.coverageoverSampdirect = coverage(thisTruth, allresoverSamp$upperdirect, allresoverSamp$lowerdirect, logit=useLogit)
    
    my.biasoverSampnaive = bias(thisTruth, allresoverSamp$logit.est, logit=useLogit, my.var=allresoverSamp$var.est)
    my.mseoverSampnaive = mse(thisTruth, allresoverSamp$logit.est, logit=useLogit, my.var=allresoverSamp$var.est)
    my.dssoverSampnaive = dss(thisTruth, allresoverSamp$logit.est, allresoverSamp$var.est)
    my.crpsoverSampnaive = crpsNormal(thisTruth, allresoverSamp$logit.est, allresoverSamp$var.est, resultType=resultType)
    my.coverageoverSampnaive = coverage(thisTruth, allresoverSamp$upper, allresoverSamp$lower, logit=useLogit)
    
    my.biasoverSampmercer = bias(thisTruth, allresoverSamp$logit.est.mercer, logit=useLogit, my.var=allresoverSamp$var.est.mercer)
    my.mseoverSampmercer = mse(thisTruth, allresoverSamp$logit.est.mercer, logit=useLogit, my.var=allresoverSamp$var.est.mercer)
    my.dssoverSampmercer = dss(thisTruth, allresoverSamp$logit.est.mercer, allresoverSamp$var.est.mercer)
    my.crpsoverSampmercer = crpsNormal(thisTruth, allresoverSamp$logit.est.mercer, allresoverSamp$var.est.mercer, resultType=resultType)
    my.coverageoverSampmercer = coverage(thisTruth, allresoverSamp$lower.mercer, allresoverSamp$upper.mercer, logit=useLogit)
    
    my.biasoverSampbymNoUrb = bias(thisTruth, designResNoUrb$overSampDat$mean[,i], logit=useLogit, my.var=(designResNoUrb$overSampDat$stddev[,i])^2)
    my.mseoverSampbymNoUrb = mse(thisTruth, designResNoUrb$overSampDat$mean[,i], logit=useLogit, my.var=(designResNoUrb$overSampDat$stddev[,i])^2)
    my.dssoverSampbymNoUrb = dss(thisTruth, designResNoUrb$overSampDat$mean[,i], (designResNoUrb$overSampDat$stddev[,i])^2)
    my.crpsoverSampbymNoUrb = crpsNormal(thisTruth, designResNoUrb$overSampDat$mean[,i], (designResNoUrb$overSampDat$stddev[,i])^2, resultType=resultType)
    my.coverageoverSampbymNoUrb = coverage(thisTruth, designResNoUrb$overSampDat$Q10[,i],designResNoUrb$overSampDat$Q90[,i], logit=useLogit)
    
    my.biasoverSampbym = bias(thisTruth, designRes$overSampDat$mean[,i], logit=useLogit, my.var=(designRes$overSampDat$stddev[,i])^2)
    my.mseoverSampbym = mse(thisTruth, designRes$overSampDat$mean[,i], logit=useLogit, my.var=(designRes$overSampDat$stddev[,i])^2)
    my.dssoverSampbym = dss(thisTruth, designRes$overSampDat$mean[,i], (designRes$overSampDat$stddev[,i])^2)
    my.crpsoverSampbym = crpsNormal(thisTruth, designRes$overSampDat$mean[,i], (designRes$overSampDat$stddev[,i])^2, resultType=resultType)
    my.coverageoverSampbym = coverage(thisTruth, designRes$overSampDat$Q10[,i],designRes$overSampDat$Q90[,i], logit=useLogit)
    
    my.biasoverSampspdeNoUrb = bias(thisTruth, allresoverSamp$logit.est.spdeNoUrb, logit=useLogit, my.var=allresoverSamp$var.est.spdeNoUrb)
    my.mseoverSampspdeNoUrb = mse(thisTruth, allresoverSamp$logit.est.spdeNoUrb, logit=useLogit, my.var=allresoverSamp$var.est.spdeNoUrb)
    my.dssoverSampspdeNoUrb = dss(thisTruth, allresoverSamp$logit.est.spdeNoUrb, allresoverSamp$var.est.spdeNoUrb)
    my.crpsoverSampspdeNoUrb = crpsNormal(thisTruth, allresoverSamp$logit.est.spdeNoUrb, allresoverSamp$var.est.spdeNoUrb, resultType=resultType)
    my.coverageoverSampspdeNoUrb = coverage(thisTruth, allresoverSamp$lower.spdeNoUrb, allresoverSamp$upper.spdeNoUrb, logit=useLogit)
    
    my.biasoverSampspde = bias(thisTruth, allresoverSamp$logit.est.spde, logit=useLogit, my.var=allresoverSamp$var.est.spde)
    my.mseoverSampspde = mse(thisTruth, allresoverSamp$logit.est.spde, logit=useLogit, my.var=allresoverSamp$var.est.spde)
    my.dssoverSampspde = dss(thisTruth, allresoverSamp$logit.est.spde, allresoverSamp$var.est.spde)
    my.crpsoverSampspde = crpsNormal(thisTruth, allresoverSamp$logit.est.spde, allresoverSamp$var.est.spde, resultType=resultType)
    my.coverageoverSampspde = coverage(thisTruth, allresoverSamp$lower.spde, allresoverSamp$upper.spde, logit=useLogit)
    
    scoresDirectoverSamp <- rbind(scoresDirectoverSamp, 
                                  data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                             bias=my.biasoverSampdirect, 
                                             mse=my.mseoverSampdirect,
                                             dss=my.dssoverSampdirect,
                                             coverage=my.coverageoverSampdirect, 
                                             var=mean(allresoverSamp$var.estdirect),
                                             crps=my.crpsoverSampdirect))
    scoresNaiveoverSamp<- rbind(scoresNaiveoverSamp, 
                                data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                           bias=my.biasoverSampnaive, 
                                           mse=my.mseoverSampnaive,
                                           dss=my.dssoverSampnaive,
                                           coverage=my.coverageoverSampnaive, 
                                           var=mean(allresoverSamp$var.est),
                                           crps=my.crpsoverSampnaive))
    scoresMerceroverSamp<- rbind(scoresMerceroverSamp, 
                                 data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                            bias=my.biasoverSampmercer, 
                                            mse=my.mseoverSampmercer,
                                            dss=my.dssoverSampmercer,
                                            coverage=my.coverageoverSampmercer, 
                                            var=mean(allresoverSamp$var.est.mercer),
                                            crps=my.crpsoverSampmercer))
    
    scoresBYMNoUrboverSamp <- rbind(scoresBYMNoUrboverSamp, 
                                    data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                               bias=my.biasoverSampbymNoUrb, 
                                               mse=my.mseoverSampbymNoUrb,
                                               dss=my.dssoverSampbymNoUrb,
                                               coverage=my.coverageoverSampbymNoUrb, 
                                               var=mean((designResNoUrb$overSampDat$stddev[,i])^2),
                                               crps=my.crpsoverSampbymNoUrb))
    
    scoresBYMoverSamp <- rbind(scoresBYMoverSamp, 
                               data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                          bias=my.biasoverSampbym, 
                                          mse=my.mseoverSampbym,
                                          dss=my.dssoverSampbym,
                                          coverage=my.coverageoverSampbym, 
                                          var=mean((designRes$overSampDat$stddev[,i])^2),
                                          crps=my.crpsoverSampbym))
    
    scoresSPDENoUrboverSamp <- rbind(scoresSPDENoUrboverSamp, 
                                     data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                                bias=my.biasoverSampspdeNoUrb, 
                                                mse=my.mseoverSampspdeNoUrb,
                                                dss=my.dssoverSampspdeNoUrb,
                                                coverage=my.coverageoverSampspdeNoUrb, 
                                                var=mean(allresoverSamp$var.est.spdeNoUrb),
                                                crps=my.crpsoverSampspdeNoUrb))
    
    scoresSPDEoverSamp <- rbind(scoresSPDEoverSamp, 
                                data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                           bias=my.biasoverSampspde, 
                                           mse=my.mseoverSampspde,
                                           dss=my.dssoverSampspde,
                                           coverage=my.coverageoverSampspde, 
                                           var=mean(allresoverSamp$var.est.spde),
                                           crps=my.crpsoverSampspde))
  }
}

if((resultType == "EA" || resultType == "pixel") && !useLogit) {
  # convert all variances to binomial scale
  N = 25
  if(sampling == "SRS") {
    scoresDirectSRS$var = scoresDirectSRS$var * N
    scoresNaiveSRS$var = scoresNaiveSRS$var * N
    scoresMercerSRS$var = scoresMercerSRS$var * N
    scoresBYMNoUrbSRS$var = scoresBYMNoUrbSRS$var * N
    scoresBYMSRS$var = scoresBYMSRS$var * N
    scoresSPDENoUrbSRS$var = scoresSPDENoUrbSRS$var * N
    scoresSPDESRS$var = scoresSPDESRS$var * N
  } else if(sampling == "oversamp") {
    scoresDirectoverSamp$var = scoresDirectoverSamp$var * N
    scoresNaiveoverSamp$var = scoresNaiveoverSamp$var * N
    scoresMerceroverSamp$var = scoresMerceroverSamp$var * N
    scoresBYMNoUrboverSamp$var = scoresBYMNoUrboverSamp$var * N
    scoresBYMoverSamp$var = scoresBYMoverSamp$var * N
    scoresSPDENoUrboverSamp$var = scoresSPDENoUrboverSamp$var * N
    scoresSPDEoverSamp$var = scoresSPDEoverSamp$var * N
  }
}

# final table (SRS):
if(sampling == "SRS") {
  naive = apply(scoresNaiveSRS[, c("bias", "mse", "dss", "crps", "var", "coverage", "length")], 2, mean)
  direct = apply(scoresDirectSRS[, c("bias", "mse", "dss", "crps", "var","coverage", "length")], 2, mean)
  mercer = apply(scoresMercerSRS[, c("bias", "mse", "dss", "crps", "var","coverage", "length")], 2, mean)
  bymNoUrb = apply(scoresBYMNoUrbSRS[, c("bias", "mse", "dss", "crps","var", "coverage", "length")], 2, mean)
  bym = apply(scoresBYMSRS[, c("bias", "mse", "dss", "crps","var", "coverage", "length")], 2, mean)
  spdeNoUrb = apply(scoresSPDENoUrbSRS[, c("bias", "mse", "dss", "crps","var", "coverage", "length")], 2, mean)
  spde = apply(scoresSPDESRS[, c("bias", "mse", "dss", "crps","var", "coverage", "length")], 2, mean)
  idx = c(1,2,5,4,6,7)
  tab = rbind(c( naive[idx]),
              c( direct[idx]),
              c( mercer[idx]),
              c( bym[idx]), 
              c( bymNoUrb[idx]), 
              c( spde[idx]), 
              c( spdeNoUrb[idx]))
  rownames(tab) = c("Naive", "Direct estimates", "Mercer et al.", "Model-based BYM", "Model-based BYM (no urban effect)", "SPDE", "SPDE (no urban effect)")
  library(xtable)
  print(xtable(tab, digits=3), 
        only.contents=TRUE, 
        include.colnames=TRUE,
        hline.after=NULL)
} else if(sampling == "oversamp") {
  # final table (oversampled):
  naive = apply(scoresNaiveoverSamp[, c("bias", "mse", "dss", "crps", "var", "coverage")], 2, mean)
  direct = apply(scoresDirectoverSamp[, c("bias", "mse", "dss", "crps", "var","coverage")], 2, mean)
  mercer = apply(scoresMerceroverSamp[, c("bias", "mse", "dss", "crps", "var","coverage")], 2, mean)
  bymNoUrb = apply(scoresBYMNoUrboverSamp[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean)
  bym = apply(scoresBYMoverSamp[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean)
  spdeNoUrb = apply(scoresSPDENoUrboverSamp[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean)
  spde = apply(scoresSPDEoverSamp[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean)
  idx = c(1,2,5,4,6)
  tab = rbind(c( naive[idx]),
              c( direct[idx]),
              c( mercer[idx]),
              c( bym[idx]), 
              c( bymNoUrb[idx]), 
              c( spde[idx]), 
              c( spdeNoUrb[idx]))
  rownames(tab) = c("Naive", "Direct estimates", "Mercer et al.", "Model-based BYM", "Model-based BYM (no urban effect)", "SPDE", "SPDE (no urban effect)")
  library(xtable)
  print(xtable(tab, digits=3), 
        only.contents=TRUE, 
        include.colnames=TRUE,
        hline.after=NULL)
}





naive = round(apply(scoresNaiveoverSamp[, c("bias","mse", "dss",  "crps", "var","coverage")], 2, mean),3)
direct = round(apply(scoresDirectoverSamp[, c("bias","mse", "dss",  "crps","var","coverage")], 2, mean),3)
mercer = round(apply(scoresMerceroverSamp[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean),3)
bym = round(apply(scoresBYMoverSamp[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean),3)
spde = round(apply(scoresSPDEoverSamp[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean),3)
tab = rbind(c( naive[idx]),
            c( direct[idx]),
            c( mercer[idx]),
            c( bym[idx]), 
            c( spde[idx]))
rownames(tab) = c("Naive", "Direct estimates", "Mercer et al.", "Model-based BYM", "SPDE")
library(xtable)
print(xtable(tab, digits=3), only.contents=TRUE, 
      include.rownames=FALSE,
      include.colnames=FALSE,
      hline.after=NULL)




naive = apply(scoresNaiveSRS[, c("bias", "mse",  "coverage")], 2, mean)
direct = apply(scoresDirectSRS[, c("bias", "mse", "coverage")], 2, mean)
mercer = apply(scoresMercerSRS[, c("bias", "mse", "coverage")], 2, mean)
bym = apply(scoresBYMSRS[, c("bias", "mse", "coverage")], 2, mean)
spde = apply(scoresSPDESRS[, c("bias", "mse", "coverage")], 2, mean)
idx = c(1,2)
tab = rbind(c( naive[idx]),
            c( direct[idx]),
            c( mercer[idx]),
            c( bym[idx]), 
            c( spde[idx]))
rownames(tab) = c("Naive", "Direct estimates", "Mercer et al.", "Model-based BYM", "SPDE")
library(xtable)
print(xtable(tab, digits=3), 
      only.contents=TRUE, 
      include.colnames=TRUE,
      hline.after=NULL)

naive = round(apply(scoresNaiveoverSamp[, c("bias","mse", "coverage")], 2, mean),3)
direct = round(apply(scoresDirectoverSamp[, c("bias","mse", "coverage")], 2, mean),3)
mercer = round(apply(scoresMerceroverSamp[, c("bias", "mse", "coverage")], 2, mean),3)
bym = round(apply(scoresBYMoverSamp[, c("bias", "mse", "coverage")], 2, mean),3)
spde = round(apply(scoresSPDEoverSamp[, c("bias", "mse", "coverage")], 2, mean),3)

tab = rbind(c( naive[idx]),
            c( direct[idx]),
            c( mercer[idx]),
            c( bym[idx]), 
            c( spde[idx]))
rownames(tab) = c("Naive", "Direct estimates", "Mercer et al.", "Model-based BYM", "SPDE")
library(xtable)
print(xtable(tab, digits=3), only.contents=TRUE, 
      include.rownames=TRUE,
      include.colnames=TRUE,
      hline.after=NULL)


# compare all four 
# pdf("Figures/biasbyregionoverSamp.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(bias~region, data=scoresDirectoverSamp, at=seq(-1, 229, by=5),
#         col="yellow", xlim=c(0,232), names=FALSE, xaxt="n")
# boxplot(bias~region, data=scoresNaiveoverSamp, at=seq(0, 230, by=5), col="orange", xlim=c(0,230), add=TRUE)
# boxplot(bias~region, data=scoresMerceroverSamp, at=seq(1, 231, by=5), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
# boxplot(bias~region, data=scoresBYMoverSamp, at=seq(2, 232, by=5), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
# axis(2, at=seq(0.5, 230.5, by=5), labels=scoresDirectoverSamp$region[1:47])
# legend("top", c("Direct estimates", "Naive", "Mercer", "BYM"),
#        fill = c("yellow", "orange", "green", "lightblue"), ncol=4, cex=2)
# abline(h=0, lwd=2, col=2)
# dev.off()
# 
# pdf("Figures/biasbyregionSRS.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(bias~region, data=scoresDirectoverSamp, at=seq(-1, 229, by=5), 
#         col="yellow", xlim=c(0,232), names=FALSE, xaxt="n")
# boxplot(bias~region, data=scoresNaiveSRS, at=seq(0, 230, by=5), col="orange", xlim=c(0,230), add=TRUE)
# boxplot(bias~region, data=scoresMercerSRS, at=seq(1, 231, by=5), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
# boxplot(bias~region, data=scoresBYMSRS, at=seq(2, 232, by=5), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
# axis(2, at=seq(0.5, 230.5, by=5), labels=scoresDirectSRS$region[1:47])
# legend("top", c("Direct estimates", "Naive", "Mercer", "BYM"),
#        fill = c("yellow", "orange", "green", "lightblue"), ncol=4, cex=2)
# abline(h=0, lwd=2, col=2)
# dev.off()
# 
# pdf("Figures/crpsbyregionoverSamp.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(crps~region, data=scoresDirectoverSamp, at=seq(-1, 229, by=5), 
#         col="yellow", xlim=c(0,232), names=FALSE, xaxt="n")
# boxplot(crps~region, data=scoresNaiveoverSamp, at=seq(0, 230, by=5), col="orange", xlim=c(0,230), add=TRUE)
# boxplot(crps~region, data=scoresMerceroverSamp, at=seq(1, 231, by=5), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
# boxplot(crps~region, data=scoresBYMoverSamp, at=seq(2, 232, by=5), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
# axis(2, at=seq(0.5, 230.5, by=5), labels=scoresDirectoversamp$region[1:47])
# legend("top", c("Direct estimates", "Naive", "Mercer", "BYM"),
#        fill = c("yellow", "orange", "green", "lightblue"), ncol=4, cex=2)
# dev.off()
# 
# pdf("Figures/crpsbyregionSRS.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(crps~region, data=scoresDirectoverSamp, at=seq(-1, 229, by=5), 
#         col="yellow", xlim=c(0,232), names=FALSE, xaxt="n")
# boxplot(crps~region, data=scoresNaiveSRS, at=seq(0, 230, by=5), col="orange", xlim=c(0,230), add=TRUE)
# boxplot(crps~region, data=scoresMercerSRS, at=seq(1, 231, by=5), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
# boxplot(crps~region, data=scoresBYMSRS, at=seq(2, 232, by=5), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
# axis(2, at=seq(0.5, 230.5, by=5), labels=scoresDirectSRS$region[1:47])
# legend("top", c("Direct estimates", "Naive", "Mercer", "BYM"),
#        fill = c("yellow", "orange", "green", "lightblue"), ncol=4, cex=2)
# dev.off()

# compare all five
pdf("figures/biasbyregionoverSamp.pdf", width=20, height=12)
par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
boxplot(bias~region, data=scoresDirectoverSamp, at=seq(-1, 275, by=6), 
        col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
boxplot(bias~region, data=scoresNaiveoverSamp, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
boxplot(bias~region, data=scoresMerceroverSamp, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
boxplot(bias~region, data=scoresBYMoverSamp, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
boxplot(bias~region, data=scoresSPDEoverSamp, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
# axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirectoverSamp$region[1:47])
# axis(1, at=seq(0.5, 276.5, by=6), labels=scoresDirectoverSamp$region[1:47])
legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
       fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
abline(h=0, lwd=2, col=2)
dev.off()

pdf("figures/biasbyregionSRS.pdf", width=20, height=12)
par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
boxplot(bias~region, data=scoresDirectoverSamp, at=seq(-1, 275, by=6), 
        col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
boxplot(bias~region, data=scoresNaiveSRS, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
boxplot(bias~region, data=scoresMercerSRS, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
boxplot(bias~region, data=scoresBYMSRS, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
boxplot(bias~region, data=scoresSPDESRS, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
# axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirectSRS$region[1:47])
# axis(1, at=seq(0.5, 276.5, by=6), labels=scoresDirectSRS$region[1:47])
legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
       fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
abline(h=0, lwd=2, col=2)
dev.off()

pdf("figures/crpsbyregionoverSamp.pdf", width=20, height=12)
par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
boxplot(crps~region, data=scoresDirectoverSamp, at=seq(-1, 275, by=6), 
        col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
boxplot(crps~region, data=scoresNaiveoverSamp, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
boxplot(crps~region, data=scoresMerceroverSamp, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
boxplot(crps~region, data=scoresBYMoverSamp, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
boxplot(crps~region, data=scoresSPDEoverSamp, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
# axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirectoverSamp$region[1:47])
legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
       fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
dev.off()

pdf("figures/crpsbyregionSRS.pdf", width=20, height=12)
par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
boxplot(crps~region, data=scoresDirectoverSamp, at=seq(-1, 275, by=6), 
        col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
boxplot(crps~region, data=scoresNaiveSRS, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
boxplot(crps~region, data=scoresMercerSRS, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
boxplot(crps~region, data=scoresBYMSRS, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
boxplot(crps~region, data=scoresSPDESRS, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
# axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirectSRS$region[1:47])
legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
       fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
dev.off()

# 
# 
# # compare naive, direct and mercer
# pdf("Figures/biasbyregionoverSamp.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(bias~region, data=scoresDirectoverSamp, at=seq(-1, 183, by=4), 
#         col="yellow", xlim=c(0,188), names=FALSE, xaxt="n")
# boxplot(bias~region, data=scoresNaiveoverSamp, at=seq(0, 184, by=4), col="orange", xlim=c(0,138), add=TRUE)
# boxplot(bias~region, data=scoresMerceroverSamp, at=seq(1, 185, by=4), col="green", xlim=c(0,138), add=TRUE, xaxt="n")
# axis(2, at=seq(0, 184, by=4), labels=scoresDirectoversamp$region[1:47])
# legend("top", c("Direct estimates", "Naive", "Mercer"),
#        fill = c("yellow", "orange", "green"), ncol=3, cex=2)
# abline(h=0, lwd=2, col=2)
# dev.off()
# 
# 
# pdf("Figures/biasbyregionSRS.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(bias~region, data=scoresDirectSRS, at=seq(-1, 183, by=4), 
#         col="yellow", xlim=c(0,188), names=FALSE, xaxt="n")
# boxplot(bias~region, data=scoresNaiveSRS, at=seq(0, 184, by=4), col="orange", xlim=c(0,138), add=TRUE)
# boxplot(bias~region, data=scoresMercerSRS, at=seq(1, 185, by=4), col="green", xlim=c(0,138), add=TRUE, xaxt="n")
# axis(2, at=seq(0, 184, by=4), labels=scoresDirectoversamp$region[1:47])
# legend("top", c("Direct estimates", "Naive", "Mercer"),
#        fill = c("yellow", "orange", "green"), ncol=3, cex=2)
# abline(h=0, lwd=2, col=2)
# dev.off()
# 
# 
# pdf("Figures/crpsbyregionoverSamp.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(crps~region, data=scoresDirectoverSamp, at=seq(-1, 183, by=4), 
#         col="yellow", xlim=c(0,188), names=FALSE, xaxt="n")
# boxplot(crps~region, data=scoresNaiveoverSamp, at=seq(0, 184, by=4), col="orange", xlim=c(0,138), add=TRUE)
# boxplot(crps~region, data=scoresMerceroverSamp, at=seq(1, 185, by=4), col="green", xlim=c(0,138), add=TRUE, xaxt="n")
# axis(2, at=seq(0, 184, by=4), labels=scoresDirectoversamp$region[1:47])
# legend("top", c("Direct estimates", "Naive", "Mercer"),
#        fill = c("yellow", "orange", "green"), ncol=3, cex=2)
# dev.off()
# 
# 
# pdf("Figures/crpsbyregionSRS.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(crps~region, data=scoresDirectSRS, at=seq(-1, 183, by=4), 
#         col="yellow", xlim=c(0,188), names=FALSE, xaxt="n")
# boxplot(crps~region, data=scoresNaiveSRS, at=seq(0, 184, by=4), col="orange", xlim=c(0,138), add=TRUE)
# boxplot(crps~region, data=scoresMercerSRS, at=seq(1, 185, by=4), col="green", xlim=c(0,138), add=TRUE, xaxt="n")
# axis(2, at=seq(0, 184, by=4), labels=scoresDirectoversamp$region[1:47])
# legend("top", c("Direct estimates", "Naive", "Mercer"),
#        fill = c("yellow", "orange", "green"), ncol=3, cex=2)
# dev.off()

# comparing only naive to direct
# 
# pdf("Figures/biasbyregionoverSamp.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(bias~region, data=scoresDirectoverSamp, at=seq(-1, 137, by=3), 
#         col="yellow", xlim=c(0,138), names=FALSE, xaxt="n")
# boxplot(bias~region, data=scoresNaiveoverSamp, at=seq(0, 138, by=3), col="orange", xlim=c(0,138), add=TRUE)
# axis(2, at=seq(-0.5, 137.5, by=3), labels=scoresDirectoversamp$region[1:47])
# legend("top", c("Direct estimates", "Naive"),
#        fill = c("yellow", "orange"), ncol=2, cex=2)
# abline(h=0, lwd=2, col=2)
# dev.off()

# pdf("Figures/biasbyregionSRS.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(bias~region, data=scoresDirectSRS, at=seq(-1, 137, by=3), 
#         col="yellow", xlim=c(0,138), names=FALSE, xaxt="n")
# boxplot(bias~region, data=scoresNaiveSRS, at=seq(0, 138, by=3), col="orange", xlim=c(0,138), add=TRUE)
# axis(2, at=seq(-0.5, 137.5, by=3), labels=scoresDirectSRS$region[1:47])
# legend("top", c("Direct estimates", "Naive"),
#        fill = c("yellow", "orange"), ncol=2, cex=2)
# abline(h=0, lwd=2, col=2)
# dev.off()

# 
# pdf("Figures/crpsbyregionoverSamp.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(crps~region, data=scoresDirectoverSamp, at=seq(-1, 137, by=3), 
#         col="yellow", xlim=c(0,138), names=FALSE, xaxt="n")
# boxplot(crps~region, data=scoresNaiveoverSamp, at=seq(0, 138, by=3), col="orange", xlim=c(0,138), add=TRUE)
# axis(2, at=seq(-0.5, 137.5, by=3), labels=scoresDirectoversamp$region[1:47])
# legend("top", c("Direct estimates", "Naive"),
#        fill = c("yellow", "orange"), ncol=2, cex=2)
# dev.off()
# 
# 
# pdf("Figures/crpsbyregionSRS.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(crps~region, data=scoresDirectSRS, at=seq(-1, 137, by=3), 
#         col="yellow", xlim=c(0,138), names=FALSE, xaxt="n")
# boxplot(crps~region, data=scoresNaiveSRS, at=seq(0, 138, by=3), col="orange", xlim=c(0,138), add=TRUE)
# axis(2, at=seq(-0.5, 137.5, by=3), labels=scoresDirectSRS$region[1:47])
# legend("top", c("Direct estimates", "Naive"),
#        fill = c("yellow", "orange"), ncol=2, cex=2)
# dev.off()