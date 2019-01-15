library(xtable)
source("scores.R")

# This model produces an xtable summarizing the scoring rules aggregated at the given level of interest
# test: whether or not to use the test data set based results or the true
#       simulation study results. The test data set results are based on fewer
#       observations, requiring fewer computations and should be used to test the code
# tausq: the spatial nugget in the simulated data, can only be the default or 0
# resultType: the level of aggregation over which scoring rules should be computed
# sampling: whether or not to use the simple random sample or urban over sampled results
# recomputeTruth: by default the should be set to TRUE, unless the user calls this function 
#                 twice in a row with the same set of inputs
# modelsI: which model results should be included. Models come in the order:
#          c("naive", "direct", "mercer", "bym", "bymNoUrb", "spde", "spdeNoUrb"), 
#          so the default value of this argument is to include results for the naive, 
#          direct, and mercer models. 1:7 is all models
# produceFigures: whether or not to produce figures based on the results
# uselogit: whether or not to produce scoring rule results at the logit or count scale. 
#           By default this is to use loaded scale results unless the resultType is EA or pixel
runCompareModels = function(test=FALSE, tausq=.1^2, resultType=c("county", "pixel", "EA"), 
                            sampling=c("SRS", "oversamp"), recomputeTruth=TRUE, modelsI=1:3, 
                            produceFigures=FALSE, uselogit=NULL) {
  # match the arguments with their correct values
  resultType = match.arg(resultType)
  sampling = match.arg(sampling)
  
  if(tausq == .1^2) {
    # out = load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOver2.RData")
    if(test) {
      out = load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOverSamplefrac0.25Test.RData")
    } else {
      out = load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOverSamplefrac0.25.RData")
    }
  }
  else {
    if(tausq != 0)
      stop("tausq can only be equal to .1^2 or 0")
    
    # out = load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOver2.RData")
    if(test) {
      out = load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOverSamplefrac0.25Test.RData")
    } else {
      out = load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOverSamplefrac0.25.RData")
    }
  }
  eaDat = SRSDat$eaDat
  
  if(sampling == "SRS") {
    clustDat = SRSDat
  } else {
    clustDat = overSampDat
  }
  
  ## generate the true values at the county level for the 
  ## given settings if these are new settings
  if(recomputeTruth){
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
  } else {
    # load the truth
    load(paste0("truthby", resultType, ".RData"))
  }
  
  ## Now pick which models to include in the results
  allModels = c("naive", "direct", "mercer", "bym", "bymNoUrb", "spde", "spdeNoUrb")
  models = allModels[modelsI]
  
  # load data
  if(tausq == .1^2) {
    if(test) {
      if("naive" %in% models || "direct" %in% models)
        out = load("resultsDirectNaiveTausq0.01Test.RData")
      if("mercer" %in% models)
        out = load("resultsMercerTausq0.01test.RData")
      if("bymNoUrb" %in% models) {
        out = load("KenyaSpatialDesignResultNewTausq0.01UrbRurFALSETest.RData")
        designResNoUrb = designRes
      }
      if("bym" %in% models)
        out = load("KenyaSpatialDesignResultNewTausq0.01UrbRurTRUETest.RData")
      if("spdeNoUrb" %in% models) {
        out = load("resultsSPDETausq0.01urbanEffectFALSETest.RData")
        spdeSRSNoUrb = spdeSRS
        spdeOverSampNoUrb = spdeOverSamp
      }
      if("spde" %in% models)
        out = load("resultsSPDETausq0.01urbanEffectTRUETest.RData")
    } else {
      if("naive" %in% models || "direct" %in% models)
        out = load("resultsDirectNaiveTausq0.01.RData")
      if("mercer" %in% models)
        out = load("resultsMercerTausq0.01.RData")
      if("bymNoUrb" %in% models) {
        out = load("KenyaSpatialDesignResultNewTausq0.01UrbRurFALSE.RData")
        designResNoUrb = designRes
      }
      if("bym" %in% models)
        out = load("KenyaSpatialDesignResultNewTausq0.01UrbRurTRUE.RData")
      if("spdeNoUrb" %in% models) {
        out = load("resultsSPDETausq0.01urbanEffectFALSE.RData")
        spdeSRSNoUrb = spdeSRS
        spdeOverSampNoUrb = spdeOverSamp
      }
      if("spde" %in% models)
        out = load("resultsSPDETausq0.01urbanEffectTRUE.RData")
    }
  } else if(tausq == 0) {
    if(test) {
      if("naive" %in% models || "direct" %in% models)
        out = load("resultsDirectNaiveTausq0Test.RData")
      if("mercer" %in% models)
        out = load("resultsMercerTausq0test.RData")
      if("bymNoUrb" %in% models) {
        out = load("KenyaSpatialDesignResultNewTausq0UrbRurFATestLSE.RData")
        designResNoUrb = designRes
      }
      if("bym" %in% models)
        out = load("KenyaSpatialDesignResultNewTausq0UrbRurTRUETest.RData")
      if("spdeNoUrb" %in% models) {
        out = load("resultsSPDETausq0urbanEffectFALSETest.RData")
        spdeSRSNoUrb = spdeSRS
        spdeOverSampNoUrb = spdeOverSamp
      }
      if("spde" %in% models)
        out = load("resultsSPDETausq0urbanEffectTRUETest.RData")
    } else {
      if("naive" %in% models || "direct" %in% models)
        out = load("resultsDirectNaiveTausq0.RData")
      if("mercer" %in% models)
        out = load("resultsMercerTausq0.RData")
      if("bymNoUrb" %in% models) {
        out = load("KenyaSpatialDesignResultNewTausq0UrbRurFALSE.RData")
        designResNoUrb = designRes
      }
      if("bym" %in% models)
        out = load("KenyaSpatialDesignResultNewTausq0UrbRurTRUE.RData")
      if("spdeNoUrb" %in% models) {
        out = load("resultsSPDETausq0urbanEffectFALSE.RData")
        spdeSRSNoUrb = spdeSRS
        spdeOverSampNoUrb = spdeOverSamp
      }
      if("spde" %in% models)
        out = load("resultsSPDETausq0urbanEffectTRUE.RData")
    }
  }
  
  # commented out the below code, since it's not clear within strata predictions are necessary
  # # modify discrete estimates stratified by urban/rural if resultType is at finer scale than county
  # # (currently, further resolving discrete estimates is only supported for the direct estimates)
  # if(resultType == "EA") {
  #   filterByUrban = function(i) {
  #     urbanI = 1:nrow(directEstSRS[[i]])
  #   }
  #   directEstSRS = lapply(1:100, )
  # }
  
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
  if("bymNoUrb" %in% models) {
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
  }
  
  # including urban effect
  if("bym" %in% models) {
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
  }
  
  for(i in c(1:100)) { # for problem fitting mercerSRS for SRS sampling, tausq=0
    # for(i in 1:100) {
    print(i)
    resultName = paste0(resultType, "Results")
    if(resultType == "EA")
      resultName = "eaResults"
    
    # convert results to the desired aggregation level
    if(sampling == "SRS") {
      if("direct" %in% models)
        directEstSRSi = getSubLevelResults(directEstSRS[[i]])
      if("naive" %in% models)
        naiveSRSi = getSubLevelResults(naiveSRS[[i]])
      if("mercer" %in% models)
        mercerSRSi = getSubLevelResults(mercerSRS[[i]])
      if(resultType != "county") {
        if("spdeNoUrb" %in% models)
          spdeSRSNoUrbi = spdeSRSNoUrb[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
        if("spde" %in% models)
          spdeSRSi = spdeSRS[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
        if("direct" %in% models)
          directEstSRSi = getSubLevelResults(directEstSRS[[i]])
      } else {
        if("spdeNoUrb" %in% models)
          spdeSRSNoUrbi = spdeSRSNoUrb[[resultName]][[i]]
        if("spde" %in% models)
          spdeSRSi = spdeSRS[[resultName]][[i]]
      }
    } else if(sampling == "oversamp") {
      if("direct" %in% models) {
        directEstoverSampi = getSubLevelResults(directEstoverSamp[[i]])
        # directEstoverSampUrbani = getSubLevelResults(directEstoverSampUrb[[i]])
        # directEstoverSampRurali = getSubLevelResults(directEstoverSampRural[[i]])
      }
      if("naive" %in% models)
        naiveoverSampi = getSubLevelResults(naiveoverSamp[[i]])
      if("mercer" %in% models)
        merceroverSampi = getSubLevelResults(merceroverSamp[[i]])
      if(resultType != "county") {
        if("spdeNoUrb" %in% models)
          spdeOverSampNoUrbi = spdeOverSampNoUrb[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
        if("spde" %in% models)
          spdeOverSampi = spdeOverSamp[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
      } else {
        if("spdeNoUrb" %in% models)
          spdeOverSampNoUrbi = spdeOverSampNoUrb[[resultName]][[i]]
        if("spde" %in% models)
          spdeOverSampi = spdeOverSamp[[resultName]][[i]]
      }
    }
    
    if(resultType == "EA") {
      # set first row of spde results to be the EA index
      if(sampling == "SRS") {
        if("spdeNoUrb" %in% models) {
          spdeSRSNoUrbi[[resultType]] = 1:nrow(spdeSRSNoUrbi)
          
          whichName = which(names(spdeSRSNoUrbi) == "EA")
          spdeSRSNoUrbi = cbind(spdeSRSNoUrbi[,whichName], spdeSRSNoUrbi[,-whichName])
          
          names(spdeSRSNoUrbi)[1] = "EA"
        }
        if("spde" %in% models) {
          spdeSRSi[[resultType]] = 1:nrow(spdeSRSi)
          
          whichName = which(names(spdeSRSi) == "EA")
          spdeSRSi = cbind(spdeSRSi[,whichName], spdeSRSi[,-whichName])
          
          names(spdeSRSi)[1] = "EA"
        }
      } else {
        if("spdeNoUrb" %in% models) {
          spdeOverSampNoUrbi[[resultType]] = 1:nrow(spdeOverSampNoUrbi)
          
          whichName = which(names(spdeOverSampNoUrbi) == "EA")
          spdeOverSampNoUrbi = cbind(spdeOverSampNoUrbi[,whichName], spdeOverSampNoUrbi[,-whichName])
          
          names(spdeOverSampNoUrbi)[1] = "EA"
        }
        if("spde" %in% models) {
          spdeOverSampi[[resultType]] = 1:nrow(spdeOverSampi)
          
          whichName = which(names(spdeOverSampi) == "EA")
          spdeOverSampi = cbind(spdeOverSampi[,whichName], spdeOverSampi[,-whichName])
          
          names(spdeOverSampi)[1] = "EA"
        }
      }
    }
    
    # for spde results, modify the name of the results
    # modify result row and column table names according to aggregation level
    if(resultType == "county") {
      # without urban effect:
      if("spdeNoUrb" %in% models) {
        if(sampling == "SRS") {
          whichName = which(names(spdeSRSNoUrbi) == "admin1")
          names(spdeSRSNoUrbi)[whichName] = resultType
        } else if (sampling == "oversamp") {
          whichName = which(names(spdeOverSampNoUrbi) == "admin1")
          names(spdeOverSampNoUrbi)[whichName] = resultType
        }
      }
      
      if("spde" %in% models) {
        # with urban effect:
        if(sampling == "SRS") {
          whichName = which(names(spdeSRSi) == "admin1")
          names(spdeSRSi)[whichName] = resultType
        } else if(sampling == "oversamp") {
          whichName = which(names(spdeOverSampi) == "admin1")
          names(spdeOverSampi)[whichName] = resultType
        }
      }
    }
    
    if(resultType == "pixel") {
      # set first row of spde results to be the pixel index
      if(sampling == "SRS") {
        if("spdeNoUrb" %in% models) {
          spdeSRSNoUrbi[[resultType]] = truth$pixel
          
          whichName = which(names(spdeSRSNoUrbi) == "pixel")
          spdeSRSNoUrbi = cbind(spdeSRSNoUrbi[,whichName], spdeSRSNoUrbi[,-whichName])
          
          names(spdeSRSNoUrbi)[1] = "pixel"
        }
        if("spde" %in% models) {
          spdeSRSi[[resultType]] = truth$pixel
          
          whichName = which(names(spdeSRSi) == "pixel")
          spdeSRSi = cbind(spdeSRSi[,whichName], spdeSRSi[,-whichName])
          
          names(spdeSRSi)[1] = "pixel"
        }
      } else {
        if("spdeNoUrb" %in% models) {
          spdeOverSampNoUrbi[[resultType]] = truth$pixel
          
          whichName = which(names(spdeOverSampNoUrbi) == "pixel")
          spdeOverSampNoUrbi = cbind(spdeOverSampNoUrbi[,whichName], spdeOverSampNoUrbi[,-whichName])
          
          names(spdeOverSampNoUrbi)[1] = "pixel"
        }
        if("spde" %in% models) {
          spdeOverSampi[[resultType]] = truth$pixel
          
          whichName = which(names(spdeOverSampi) == "pixel")
          spdeOverSampi = cbind(spdeOverSampi[,whichName], spdeOverSampi[,-whichName])
          
          names(spdeOverSampi)[1] = "pixel"
        }
      }
    }
    
    # change names of table variables in spde model with no urban effect to reflect that
    if(sampling == "SRS") {
      if("spdeNoUrb" %in% models)
        names(spdeSRSNoUrbi)[2:6] = paste0(names(spdeSRSNoUrbi)[2:6], "NoUrb")
    } else {
      if("spdeNoUrb" %in% models)
        names(spdeOverSampNoUrbi)[2:6] = paste0(names(spdeOverSampNoUrbi)[2:6], "NoUrb")
    }
    
    if(sampling == "SRS") {
      if("direct" %in% models) {
        allresSRS = merge(truth, directEstSRSi, by=resultType)
        colnames(allresSRS) = c(resultType, "truth", paste(colnames(allresSRS)[3:8], "direct", sep=""))
      } else {
        stop("direct estimates must be included at this point in order to name the estimate table columns")
      }
      if("naive" %in% models)
        allresSRS = merge(allresSRS, naiveSRSi, by=resultType)
      if("mercer" %in% models)
        allresSRS = merge(allresSRS, mercerSRSi, by=resultType)
      if("spdeNoUrb" %in% models)
        allresSRS = merge(allresSRS, spdeSRSNoUrbi, by=resultType)
      if("spde" %in% models)
        allresSRS = merge(allresSRS, spdeSRSi, by=resultType)
    }
    
    if(sampling == "oversamp") {
      if("direct" %in% models) {
        allresoverSamp = merge(truth, directEstoverSampi, by=resultType)
        colnames(allresoverSamp) = c(resultType, "truth", paste(colnames(allresoverSamp)[3:8], "direct", sep=""))
      } else {
        stop("direct estimates must be included at this point in order to name the estimate table columns")
      }
      if("knife" %in% models)
        allresoverSamp = merge(allresoverSamp, naiveoverSampi, by=resultType)
      if("mercer" %in% models)
        allresoverSamp = merge(allresoverSamp, merceroverSampi, by=resultType)
      if("spdeNoUrb" %in% models)
        allresoverSamp = merge(allresoverSamp, spdeOverSampNoUrbi, by=resultType)
      if("spde" %in% models)
        allresoverSamp = merge(allresoverSamp, spdeOverSampi, by=resultType)
    }
    
    # set whether or not to calculate scores on logit scale depending on result type
    thisTruth = allresSRS$truth
    if(is.null(uselogit)) {
      useLogit = FALSE
      if(resultType != "EA" && resultType != "pixel") {
        thisTruth = logit(thisTruth)
        useLogit=TRUE
      }
    }
    
    # SRS setting 
    if(sampling == "SRS") {
      if("direct" %in% models) {
        my.biasSRSdirect = bias(thisTruth, allresSRS$logit.estdirect, logit=useLogit, my.var=allresSRS$var.estdirect)
        my.mseSRSdirect = mse(thisTruth, allresSRS$logit.estdirect, logit=useLogit, my.var=allresSRS$var.estdirect)
        my.dssSRSdirect = dss(thisTruth, allresSRS$logit.estdirect, allresSRS$var.estdirect)
        my.crpsSRSdirect = crpsNormal(thisTruth, allresSRS$logit.estdirect, allresSRS$var.estdirect, resultType=resultType)
        my.coverageSRSdirect = coverage(thisTruth, allresSRS$lowerdirect, allresSRS$upperdirect, logit=useLogit)
        my.lengthSRSdirect = intervalWidth(allresSRS$lowerdirect, allresSRS$upperdirect, logit=useLogit)
      }
      
      if("naive" %in% models) {
        my.biasSRSnaive = bias(thisTruth, allresSRS$logit.est, logit=useLogit, my.var=allresSRS$var.est)
        my.mseSRSnaive = mse(thisTruth, allresSRS$logit.est, logit=useLogit, my.var=allresSRS$var.est)
        my.dssSRSnaive = dss(thisTruth, allresSRS$logit.est, allresSRS$var.est)
        my.crpsSRSnaive = crpsNormal(thisTruth, allresSRS$logit.est, allresSRS$var.est, resultType=resultType)
        my.coverageSRSnaive = coverage(thisTruth, allresSRS$lower, allresSRS$upper, logit=useLogit)
        my.lengthSRSnaive = intervalWidth(allresSRS$lower, allresSRS$upper, logit=useLogit)
      }
      
      if("mercer" %in% models) {
        my.biasSRSmercer = bias(thisTruth, allresSRS$logit.est.mercer, logit=useLogit, my.var=allresSRS$var.est.mercer)
        my.mseSRSmercer = mse(thisTruth, allresSRS$logit.est.mercer, logit=useLogit, my.var=allresSRS$var.est.mercer)
        my.dssSRSmercer = dss(thisTruth, allresSRS$logit.est.mercer, allresSRS$var.est.mercer)
        my.crpsSRSmercer = crpsNormal(thisTruth, allresSRS$logit.est.mercer, allresSRS$var.est.mercer, resultType=resultType)
        my.coverageSRSmercer = coverage(thisTruth, allresSRS$lower.mercer, allresSRS$upper.mercer, logit=useLogit)
        my.lengthSRSmercer = intervalWidth(allresSRS$lower.mercer, allresSRS$upper.mercer, logit=useLogit)
      }
      
      if("bymNoUrb" %in% models) {
        my.biasSRSbymNoUrb = bias(thisTruth, designResNoUrb$SRSdat$mean[,i], logit=useLogit, my.var=(designResNoUrb$SRSdat$stddev[,i])^2)
        my.mseSRSbymNoUrb = mse(thisTruth, designResNoUrb$SRSdat$mean[,i], logit=useLogit, my.var=(designResNoUrb$SRSdat$stddev[,i])^2)
        my.dssSRSbymNoUrb = dss(thisTruth, designResNoUrb$SRSdat$mean[,i], (designResNoUrb$SRSdat$stddev[,i])^2)
        my.crpsSRSbymNoUrb = crpsNormal(thisTruth, designResNoUrb$SRSdat$mean[,i], (designResNoUrb$SRSdat$stddev[,i])^2, resultType=resultType)
        my.coverageSRSbymNoUrb = coverage(thisTruth, designResNoUrb$SRSdat$Q10[,i],designResNoUrb$SRSdat$Q90[,i], logit=useLogit)
        my.lengthSRSbymNoUrb = intervalWidth(designResNoUrb$SRSdat$Q10[,i], designResNoUrb$SRSdat$Q90[,i], logit=useLogit)
      }
      
      if("bym" %in% models) {
        my.biasSRSbym = bias(thisTruth, designRes$SRSdat$mean[,i], logit=useLogit, my.var=(designRes$SRSdat$stddev[,i])^2)
        my.mseSRSbym = mse(thisTruth, designRes$SRSdat$mean[,i], logit=useLogit, my.var=(designRes$SRSdat$stddev[,i])^2)
        my.dssSRSbym = dss(thisTruth, designRes$SRSdat$mean[,i], (designRes$SRSdat$stddev[,i])^2)
        my.crpsSRSbym = crpsNormal(thisTruth, designRes$SRSdat$mean[,i], (designRes$SRSdat$stddev[,i])^2, resultType=resultType)
        my.coverageSRSbym = coverage(thisTruth, designRes$SRSdat$Q10[,i],designRes$SRSdat$Q90[,i], logit=useLogit)
        my.lengthSRSbym = intervalWidth(designRes$SRSdat$Q10[,i], designRes$SRSdat$Q90[,i], logit=useLogit)
      }
      
      if("spdeNoUrb" %in% models) {
        my.biasSRSspdeNoUrb = bias(thisTruth, allresSRS$logit.est.spdeNoUrb, logit=useLogit, my.var=allresSRS$var.est.spdeNoUrb)
        my.mseSRSspdeNoUrb = mse(thisTruth, allresSRS$logit.est.spdeNoUrb, logit=useLogit, my.var=allresSRS$var.est.spdeNoUrb)
        my.dssSRSspdeNoUrb = dss(thisTruth, allresSRS$logit.est.spdeNoUrb, allresSRS$var.est.spdeNoUrb)
        # my.crpsSRSspdeNoUrb = crpsNormal(thisTruth, allresSRS$logit.est.spdeNoUrb, allresSRS$var.est.spdeNoUrb, resultType=resultType)
        my.crpsSRSspdeNoUrb = allresSRS$crps.spdeNoUrb
        my.coverageSRSspdeNoUrb = coverage(thisTruth, allresSRS$lower.spdeNoUrb, allresSRS$upper.spdeNoUrb, logit=useLogit)
        my.lengthSRSspdeNoUrb = intervalWidth(allresSRS$lower.spdeNoUrb, allresSRS$upper.spdeNoUrb, logit=useLogit)
      }
      
      if("spde" %in% models) {
        my.biasSRSspde = bias(thisTruth, allresSRS$logit.est.spde, logit=useLogit, my.var=allresSRS$var.est.spde)
        my.mseSRSspde = mse(thisTruth, allresSRS$logit.est.spde, logit=useLogit, my.var=allresSRS$var.est.spde)
        my.dssSRSspde = dss(thisTruth, allresSRS$logit.est.spde, allresSRS$var.est.spde)
        # my.crpsSRSspde = crpsNormal(thisTruth, allresSRS$logit.est.spde, allresSRS$var.est.spde, resultType=resultType)
        my.crpsSRSspde = allresSRS$crps.spde
        my.coverageSRSspde = coverage(thisTruth, allresSRS$lower.spde, allresSRS$upper.spde, logit=useLogit)
        my.lengthSRSspde = intervalWidth(allresSRS$lower.spde, allresSRS$upper.spde, logit=useLogit)
      }
      
      if("direct" %in% models) {
        scoresDirectSRS <- rbind(scoresDirectSRS,
                                 data.frame(dataset=i, region=allresSRS[[resultType]],
                                            bias=my.biasSRSdirect,
                                            mse=my.mseSRSdirect,
                                            dss=my.dssSRSdirect,
                                            coverage=my.coverageSRSdirect,
                                            var=mean(allresSRS$var.estdirect),
                                            crps=my.crpsSRSdirect, 
                                            length=my.lengthSRSdirect))
      }
      
      if("naive" %in% models) {
        scoresNaiveSRS <- rbind(scoresNaiveSRS,
                                data.frame(dataset=i, region=allresSRS[[resultType]],
                                           bias=my.biasSRSnaive,
                                           mse=my.mseSRSnaive,
                                           dss=my.dssSRSnaive,
                                           coverage=my.coverageSRSnaive,
                                           var=mean(allresSRS$var.est),
                                           crps=my.crpsSRSnaive, 
                                           length=my.lengthSRSnaive))
      }
      
      if("mercer" %in% models) {
        scoresMercerSRS <- rbind(scoresMercerSRS,
                                 data.frame(dataset=i, region=allresSRS[[resultType]],
                                            bias=my.biasSRSmercer,
                                            mse=my.mseSRSmercer,
                                            dss=my.dssSRSmercer,
                                            coverage=my.coverageSRSmercer,
                                            var=mean(allresSRS$var.est.mercer),
                                            crps=my.crpsSRSmercer, 
                                            length=my.lengthSRSmercer))
      }
      
      if("bymNoUrb" %in% models) {
        scoresBYMNoUrbSRS <- rbind(scoresBYMNoUrbSRS,
                                   data.frame(dataset=i, region=allresSRS[[resultType]],
                                              bias=my.biasSRSbymNoUrb,
                                              mse=my.mseSRSbymNoUrb,
                                              dss=my.dssSRSbymNoUrb,
                                              coverage=my.coverageSRSbymNoUrb,
                                              var=mean((designResNoUrb$SRSdat$stddev[,i])^2),
                                              crps=my.crpsSRSbymNoUrb, 
                                              length=my.lengthSRSbymNoUrb))
      }
      
      if("bym" %in% models) {
        scoresBYMSRS <- rbind(scoresBYMSRS,
                              data.frame(dataset=i, region=allresSRS[[resultType]],
                                         bias=my.biasSRSbym,
                                         mse=my.mseSRSbym,
                                         dss=my.dssSRSbym,
                                         coverage=my.coverageSRSbym,
                                         var=mean((designRes$SRSdat$stddev[,i])^2),
                                         crps=my.crpsSRSbym, 
                                         length=my.lengthSRSbym))
      }
      
      if("spdeNoUrb" %in% models) {
        scoresSPDENoUrbSRS <- rbind(scoresSPDENoUrbSRS, 
                                    data.frame(dataset=i, region=allresSRS[[resultType]], 
                                               bias=my.biasSRSspdeNoUrb, 
                                               mse=my.mseSRSspdeNoUrb,
                                               dss=my.dssSRSspdeNoUrb,
                                               coverage=my.coverageSRSspdeNoUrb, 
                                               var=mean(allresSRS$var.est.spdeNoUrb),
                                               crps=my.crpsSRSspdeNoUrb, 
                                               length=my.lengthSRSspdeNoUrb))
      }
      
      if("spde" %in% models) {
        scoresSPDESRS <- rbind(scoresSPDESRS, 
                               data.frame(dataset=i, region=allresSRS[[resultType]], 
                                          bias=my.biasSRSspde, 
                                          mse=my.mseSRSspde,
                                          dss=my.dssSRSspde,
                                          coverage=my.coverageSRSspde, 
                                          var=mean(allresSRS$var.est.spde),
                                          crps=my.crpsSRSspde, 
                                          length=my.lengthSRSspde))
      }
      
    } else if(sampling == "oversamp") {
      # oversampling setting
      if("direct" %in% models) {
        my.biasoverSampdirect = bias(thisTruth, allresoverSamp$logit.estdirect, logit=useLogit, my.var=allresoverSamp$var.estdirect)
        my.mseoverSampdirect = mse(thisTruth, allresoverSamp$logit.estdirect, logit=useLogit, my.var=allresoverSamp$var.estdirect)
        my.dssoverSampdirect = dss(thisTruth, allresoverSamp$logit.estdirect, allresoverSamp$var.estdirect)
        my.crpsoverSampdirect = crpsNormal(thisTruth, allresoverSamp$logit.estdirect, allresoverSamp$var.estdirect, resultType=resultType)
        my.coverageoverSampdirect = coverage(thisTruth, allresoverSamp$upperdirect, allresoverSamp$lowerdirect, logit=useLogit)
      }
      
      if("naive" %in% models) {
        my.biasoverSampnaive = bias(thisTruth, allresoverSamp$logit.est, logit=useLogit, my.var=allresoverSamp$var.est)
        my.mseoverSampnaive = mse(thisTruth, allresoverSamp$logit.est, logit=useLogit, my.var=allresoverSamp$var.est)
        my.dssoverSampnaive = dss(thisTruth, allresoverSamp$logit.est, allresoverSamp$var.est)
        my.crpsoverSampnaive = crpsNormal(thisTruth, allresoverSamp$logit.est, allresoverSamp$var.est, resultType=resultType)
        my.coverageoverSampnaive = coverage(thisTruth, allresoverSamp$upper, allresoverSamp$lower, logit=useLogit)
      }
      
      if("mercer" %in% models) {
        my.biasoverSampmercer = bias(thisTruth, allresoverSamp$logit.est.mercer, logit=useLogit, my.var=allresoverSamp$var.est.mercer)
        my.mseoverSampmercer = mse(thisTruth, allresoverSamp$logit.est.mercer, logit=useLogit, my.var=allresoverSamp$var.est.mercer)
        my.dssoverSampmercer = dss(thisTruth, allresoverSamp$logit.est.mercer, allresoverSamp$var.est.mercer)
        my.crpsoverSampmercer = crpsNormal(thisTruth, allresoverSamp$logit.est.mercer, allresoverSamp$var.est.mercer, resultType=resultType)
        my.coverageoverSampmercer = coverage(thisTruth, allresoverSamp$lower.mercer, allresoverSamp$upper.mercer, logit=useLogit)
      }
      
      if("bymNoUrb" %in% models) {
        my.biasoverSampbymNoUrb = bias(thisTruth, designResNoUrb$overSampDat$mean[,i], logit=useLogit, my.var=(designResNoUrb$overSampDat$stddev[,i])^2)
        my.mseoverSampbymNoUrb = mse(thisTruth, designResNoUrb$overSampDat$mean[,i], logit=useLogit, my.var=(designResNoUrb$overSampDat$stddev[,i])^2)
        my.dssoverSampbymNoUrb = dss(thisTruth, designResNoUrb$overSampDat$mean[,i], (designResNoUrb$overSampDat$stddev[,i])^2)
        my.crpsoverSampbymNoUrb = crpsNormal(thisTruth, designResNoUrb$overSampDat$mean[,i], (designResNoUrb$overSampDat$stddev[,i])^2, resultType=resultType)
        my.coverageoverSampbymNoUrb = coverage(thisTruth, designResNoUrb$overSampDat$Q10[,i],designResNoUrb$overSampDat$Q90[,i], logit=useLogit)
      }
      
      if("bym" %in% models) {
        my.biasoverSampbym = bias(thisTruth, designRes$overSampDat$mean[,i], logit=useLogit, my.var=(designRes$overSampDat$stddev[,i])^2)
        my.mseoverSampbym = mse(thisTruth, designRes$overSampDat$mean[,i], logit=useLogit, my.var=(designRes$overSampDat$stddev[,i])^2)
        my.dssoverSampbym = dss(thisTruth, designRes$overSampDat$mean[,i], (designRes$overSampDat$stddev[,i])^2)
        my.crpsoverSampbym = crpsNormal(thisTruth, designRes$overSampDat$mean[,i], (designRes$overSampDat$stddev[,i])^2, resultType=resultType)
        my.coverageoverSampbym = coverage(thisTruth, designRes$overSampDat$Q10[,i],designRes$overSampDat$Q90[,i], logit=useLogit)
      }
      
      if("spdeNoUrb" %in% models) {
        my.biasoverSampspdeNoUrb = bias(thisTruth, allresoverSamp$logit.est.spdeNoUrb, logit=useLogit, my.var=allresoverSamp$var.est.spdeNoUrb)
        my.mseoverSampspdeNoUrb = mse(thisTruth, allresoverSamp$logit.est.spdeNoUrb, logit=useLogit, my.var=allresoverSamp$var.est.spdeNoUrb)
        my.dssoverSampspdeNoUrb = dss(thisTruth, allresoverSamp$logit.est.spdeNoUrb, allresoverSamp$var.est.spdeNoUrb)
        my.crpsoverSampspdeNoUrb = crpsNormal(thisTruth, allresoverSamp$logit.est.spdeNoUrb, allresoverSamp$var.est.spdeNoUrb, resultType=resultType)
        my.coverageoverSampspdeNoUrb = coverage(thisTruth, allresoverSamp$lower.spdeNoUrb, allresoverSamp$upper.spdeNoUrb, logit=useLogit)
      }
      
      if("spde" %in% models) {
        my.biasoverSampspde = bias(thisTruth, allresoverSamp$logit.est.spde, logit=useLogit, my.var=allresoverSamp$var.est.spde)
        my.mseoverSampspde = mse(thisTruth, allresoverSamp$logit.est.spde, logit=useLogit, my.var=allresoverSamp$var.est.spde)
        my.dssoverSampspde = dss(thisTruth, allresoverSamp$logit.est.spde, allresoverSamp$var.est.spde)
        my.crpsoverSampspde = crpsNormal(thisTruth, allresoverSamp$logit.est.spde, allresoverSamp$var.est.spde, resultType=resultType)
        my.coverageoverSampspde = coverage(thisTruth, allresoverSamp$lower.spde, allresoverSamp$upper.spde, logit=useLogit)
      }
      
      if("direct" %in% models) {
        scoresDirectoverSamp <- rbind(scoresDirectoverSamp, 
                                      data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                                 bias=my.biasoverSampdirect, 
                                                 mse=my.mseoverSampdirect,
                                                 dss=my.dssoverSampdirect,
                                                 coverage=my.coverageoverSampdirect, 
                                                 var=mean(allresoverSamp$var.estdirect),
                                                 crps=my.crpsoverSampdirect))
      }
      
      if("naive" %in% models) {
        scoresNaiveoverSamp<- rbind(scoresNaiveoverSamp, 
                                    data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                               bias=my.biasoverSampnaive, 
                                               mse=my.mseoverSampnaive,
                                               dss=my.dssoverSampnaive,
                                               coverage=my.coverageoverSampnaive, 
                                               var=mean(allresoverSamp$var.est),
                                               crps=my.crpsoverSampnaive))
      }
      
      if("mercer" %in% models) {
        scoresMerceroverSamp<- rbind(scoresMerceroverSamp, 
                                     data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                                bias=my.biasoverSampmercer, 
                                                mse=my.mseoverSampmercer,
                                                dss=my.dssoverSampmercer,
                                                coverage=my.coverageoverSampmercer, 
                                                var=mean(allresoverSamp$var.est.mercer),
                                                crps=my.crpsoverSampmercer))
      }
      
      if("bymNoUrb" %in% models) {
        scoresBYMNoUrboverSamp <- rbind(scoresBYMNoUrboverSamp, 
                                        data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                                   bias=my.biasoverSampbymNoUrb, 
                                                   mse=my.mseoverSampbymNoUrb,
                                                   dss=my.dssoverSampbymNoUrb,
                                                   coverage=my.coverageoverSampbymNoUrb, 
                                                   var=mean((designResNoUrb$overSampDat$stddev[,i])^2),
                                                   crps=my.crpsoverSampbymNoUrb))
      }
      
      if("bym" %in% models) {
        scoresBYMoverSamp <- rbind(scoresBYMoverSamp, 
                                   data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                              bias=my.biasoverSampbym, 
                                              mse=my.mseoverSampbym,
                                              dss=my.dssoverSampbym,
                                              coverage=my.coverageoverSampbym, 
                                              var=mean((designRes$overSampDat$stddev[,i])^2),
                                              crps=my.crpsoverSampbym))
      }
      
      if("spdeNoUrb" %in% models) {
        scoresSPDENoUrboverSamp <- rbind(scoresSPDENoUrboverSamp, 
                                         data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                                    bias=my.biasoverSampspdeNoUrb, 
                                                    mse=my.mseoverSampspdeNoUrb,
                                                    dss=my.dssoverSampspdeNoUrb,
                                                    coverage=my.coverageoverSampspdeNoUrb, 
                                                    var=mean(allresoverSamp$var.est.spdeNoUrb),
                                                    crps=my.crpsoverSampspdeNoUrb))
      }
      
      
      if("spde" %in% models) {
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
  }
  
  # the below code is commented out, since we want all scoring rules on the empirical proportion scale
  # if((resultType == "EA" || resultType == "pixel") && !useLogit) {
  #   # convert all variances to binomial scale
  #   N = 25
  #   if(sampling == "SRS") {
  #     scoresDirectSRS$var = scoresDirectSRS$var * N
  #     scoresNaiveSRS$var = scoresNaiveSRS$var * N
  #     scoresMercerSRS$var = scoresMercerSRS$var * N
  #     scoresBYMNoUrbSRS$var = scoresBYMNoUrbSRS$var * N
  #     scoresBYMSRS$var = scoresBYMSRS$var * N
  #     scoresSPDENoUrbSRS$var = scoresSPDENoUrbSRS$var * N
  #     scoresSPDESRS$var = scoresSPDESRS$var * N
  #   } else if(sampling == "oversamp") {
  #     scoresDirectoverSamp$var = scoresDirectoverSamp$var * N
  #     scoresNaiveoverSamp$var = scoresNaiveoverSamp$var * N
  #     scoresMerceroverSamp$var = scoresMerceroverSamp$var * N
  #     scoresBYMNoUrboverSamp$var = scoresBYMNoUrboverSamp$var * N
  #     scoresBYMoverSamp$var = scoresBYMoverSamp$var * N
  #     scoresSPDENoUrboverSamp$var = scoresSPDENoUrboverSamp$var * N
  #     scoresSPDEoverSamp$var = scoresSPDEoverSamp$var * N
  #   }
  # }
  
  # final table (SRS):
  if(sampling == "SRS") {
    if("naive" %in% models)
      naive = apply(scoresNaiveSRS[, c("bias", "mse", "dss", "crps", "var", "coverage", "length")], 2, mean)
    if("direct" %in% models)
      direct = apply(scoresDirectSRS[, c("bias", "mse", "dss", "crps", "var","coverage", "length")], 2, mean)
    if("mercer" %in% models)
      mercer = apply(scoresMercerSRS[, c("bias", "mse", "dss", "crps", "var","coverage", "length")], 2, mean)
    if("bymNoUrb" %in% models)
      bymNoUrb = apply(scoresBYMNoUrbSRS[, c("bias", "mse", "dss", "crps","var", "coverage", "length")], 2, mean)
    if("bym" %in% models)
      bym = apply(scoresBYMSRS[, c("bias", "mse", "dss", "crps","var", "coverage", "length")], 2, mean)
    if("spdeNoUrb" %in% models)
      spdeNoUrb = apply(scoresSPDENoUrbSRS[, c("bias", "mse", "dss", "crps","var", "coverage", "length")], 2, mean)
    if("spde" %in% models)
      spde = apply(scoresSPDESRS[, c("bias", "mse", "dss", "crps","var", "coverage", "length")], 2, mean)
    idx = c(1,2,5,4,6,7)
    tab = c()
    if("naive" %in% models)
      tab = rbind(tab, c(naive[idx]))
    if("direct" %in% models)
      tab = rbind(tab, c(direct[idx]))
    if("mercer" %in% models)
      tab = rbind(tab, c(mercer[idx]))
    if("bym" %in% models)
      tab = rbind(tab, c(bym[idx]))
    if("bymNoUrb" %in% models)
      tab = rbind(tab, c(bymNoUrb[idx]))
    if("spde" %in% models)
      tab = rbind(tab, c(spde[idx]))
    if("spdeNoUrb" %in% models)
      tab = rbind(tab, c(spdeNoUrb[idx]))
    allNames = c("Naive", "Direct estimates", "Mercer et al.", "Model-based BYM", "Model-based BYM (no urban effect)", "SPDE", "SPDE (no urban effect)")
    rownames(tab) = allNames[modelsI]
    
    print(xtable(tab, digits=3), 
          only.contents=TRUE, 
          include.colnames=TRUE,
          hline.after=NULL)
  } else if(sampling == "oversamp") {
    # final table (oversampled):
    if("naive" %in% models)
      naive = apply(scoresNaiveoverSamp[, c("bias", "mse", "dss", "crps", "var", "coverage", "length")], 2, mean)
    if("direct" %in% models)
      direct = apply(scoresDirectoverSamp[, c("bias", "mse", "dss", "crps", "var","coverage", "length")], 2, mean)
    if("mercer" %in% models)
      mercer = apply(scoresMerceroverSamp[, c("bias", "mse", "dss", "crps", "var","coverage", "length")], 2, mean)
    if("bymNoUrb" %in% models)
      bymNoUrb = apply(scoresBYMNoUrboverSamp[, c("bias", "mse", "dss", "crps","var", "coverage", "length")], 2, mean)
    if("bym" %in% models)
      bym = apply(scoresBYMoverSamp[, c("bias", "mse", "dss", "crps","var", "coverage", "length")], 2, mean)
    if("spdeNoUrb" %in% models)
      spdeNoUrb = apply(scoresSPDENoUrboverSamp[, c("bias", "mse", "dss", "crps","var", "coverage", "length")], 2, mean)
    if("spde" %in% models)
      spde = apply(scoresSPDEoverSamp[, c("bias", "mse", "dss", "crps","var", "coverage", "length")], 2, mean)
    idx = c(1,2,5,4,6,7)
    tab = c()
    if("direct" %in% models)
      tab = rbind(tab, c(direct[idx]))
    if("naive" %in% models)
      tab = rbind(tab, c(naive[idx]))
    if("mercer" %in% models)
      tab = rbind(tab, c(mercer[idx]))
    if("bym" %in% models)
      tab = rbind(tab, c(bym[idx]))
    if("bymNoUrb" %in% models)
      tab = rbind(tab, c(bymNoUrb[idx]))
    if("spde" %in% models)
      tab = rbind(tab, c(spde[idx]))
    if("spdeNoUrb" %in% models)
      tab = rbind(tab, c(spdeNoUrb[idx]))
    allNames = c("Naive", "Direct estimates", "Mercer et al.", "Model-based BYM", "Model-based BYM (no urban effect)", "SPDE", "SPDE (no urban effect)")
    rownames(tab) = allNames[modelsI]
    
    print(xtable(tab, digits=3), 
          only.contents=TRUE, 
          include.colnames=TRUE,
          hline.after=NULL)
    
    
    # the below code was phased out after including optional models
    # naive = apply(scoresNaiveoverSamp[, c("bias", "mse", "dss", "crps", "var", "coverage")], 2, mean)
    # direct = apply(scoresDirectoverSamp[, c("bias", "mse", "dss", "crps", "var","coverage")], 2, mean)
    # mercer = apply(scoresMerceroverSamp[, c("bias", "mse", "dss", "crps", "var","coverage")], 2, mean)
    # bymNoUrb = apply(scoresBYMNoUrboverSamp[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean)
    # bym = apply(scoresBYMoverSamp[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean)
    # spdeNoUrb = apply(scoresSPDENoUrboverSamp[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean)
    # spde = apply(scoresSPDEoverSamp[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean)
    # idx = c(1,2,5,4,6)
    # tab = rbind(c( naive[idx]),
    #             c( direct[idx]),
    #             c( mercer[idx]),
    #             c( bym[idx]), 
    #             c( bymNoUrb[idx]), 
    #             c( spde[idx]), 
    #             c( spdeNoUrb[idx]))
    # rownames(tab) = c("Naive", "Direct estimates", "Mercer et al.", "Model-based BYM", "Model-based BYM (no urban effect)", "SPDE", "SPDE (no urban effect)")
    # library(xtable)
    # print(xtable(tab, digits=3), 
    #       only.contents=TRUE, 
    #       include.colnames=TRUE,
    #       hline.after=NULL)
  }
  
  if(produceFigures) {
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
  }
  
  tab
}







##### Now run the rest of the script (except possibly for the plotting, depending 
##### on the result aggregation level).















# naive = round(apply(scoresNaiveoverSamp[, c("bias","mse", "dss",  "crps", "var","coverage")], 2, mean),3)
# direct = round(apply(scoresDirectoverSamp[, c("bias","mse", "dss",  "crps","var","coverage")], 2, mean),3)
# mercer = round(apply(scoresMerceroverSamp[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean),3)
# bym = round(apply(scoresBYMoverSamp[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean),3)
# spde = round(apply(scoresSPDEoverSamp[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean),3)
# tab = rbind(c( naive[idx]),
#             c( direct[idx]),
#             c( mercer[idx]),
#             c( bym[idx]), 
#             c( spde[idx]))
# rownames(tab) = c("Naive", "Direct estimates", "Mercer et al.", "Model-based BYM", "SPDE")
# library(xtable)
# print(xtable(tab, digits=3), only.contents=TRUE, 
#       include.rownames=FALSE,
#       include.colnames=FALSE,
#       hline.after=NULL)
# 
# 
# 
# 
# naive = apply(scoresNaiveSRS[, c("bias", "mse",  "coverage")], 2, mean)
# direct = apply(scoresDirectSRS[, c("bias", "mse", "coverage")], 2, mean)
# mercer = apply(scoresMercerSRS[, c("bias", "mse", "coverage")], 2, mean)
# bym = apply(scoresBYMSRS[, c("bias", "mse", "coverage")], 2, mean)
# spde = apply(scoresSPDESRS[, c("bias", "mse", "coverage")], 2, mean)
# idx = c(1,2)
# tab = rbind(c( naive[idx]),
#             c( direct[idx]),
#             c( mercer[idx]),
#             c( bym[idx]), 
#             c( spde[idx]))
# rownames(tab) = c("Naive", "Direct estimates", "Mercer et al.", "Model-based BYM", "SPDE")
# library(xtable)
# print(xtable(tab, digits=3), 
#       only.contents=TRUE, 
#       include.colnames=TRUE,
#       hline.after=NULL)
# 
# naive = round(apply(scoresNaiveoverSamp[, c("bias","mse", "coverage")], 2, mean),3)
# direct = round(apply(scoresDirectoverSamp[, c("bias","mse", "coverage")], 2, mean),3)
# mercer = round(apply(scoresMerceroverSamp[, c("bias", "mse", "coverage")], 2, mean),3)
# bym = round(apply(scoresBYMoverSamp[, c("bias", "mse", "coverage")], 2, mean),3)
# spde = round(apply(scoresSPDEoverSamp[, c("bias", "mse", "coverage")], 2, mean),3)
# 
# tab = rbind(c( naive[idx]),
#             c( direct[idx]),
#             c( mercer[idx]),
#             c( bym[idx]), 
#             c( spde[idx]))
# rownames(tab) = c("Naive", "Direct estimates", "Mercer et al.", "Model-based BYM", "SPDE")
# library(xtable)
# print(xtable(tab, digits=3), only.contents=TRUE, 
#       include.rownames=TRUE,
#       include.colnames=TRUE,
#       hline.after=NULL)


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