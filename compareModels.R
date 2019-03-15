library(xtable)
source("scores.R")

# This model produces an xtable summarizing the scoring rules aggregated at the given level of interest. 
# Same as runCompareModels, except uses the updated scoring rule function getScores, so it produces 
# scoring rules based on the predictive distribution with and without binomial variation.
# test: whether or not to use the test data set based results or the true
#       simulation study results. The test data set results are based on fewer
#       observations, requiring fewer computations and should be used to test the code
# tausq: the spatial nugget in the simulated data, can only be the default or 0
# resultType: the level of aggregation over which scoring rules should be computed
# sampling: whether or not to use the simple random sample or urban over sampled results
# recomputeTruth: by default the should be set to TRUE, unless the user calls this function 
#                 twice in a row with the same set of inputs
# modelsI: which model results should be included. Models come in the order:
#          c("naive", "direct", "mercer", "bym", "bymMod", "bymNoUrb", "bymNoUrbMod", "bymNoClust", "bymNoUrbClust", "spde", "spdeNoUrb"), 
#          so the default value of this argument is to include results for the naive, 
#          direct, and mercer models. 1:11 is all models
# produceFigures: whether or not to produce figures based on the results
# printIEvery: how often to print progress
# maxDataSets: if not null, set this to a small integer to only score this many datasets 
#              for testing purposes
# nsim: the number of points with which to approximate the distribution of probability on the logit scale
# saveResults: whether or not to save a disk image of the results
# loadResults: if TRUE loads the saved disk image rather than recomputing results
# xtable.args: the arguments passed to xtable for printing results
# tableFormat: If "1", binomial scores are considered the same model as non-binomial scores. 
#              If "2", binomial scores are put on extra rows of the printed table
runCompareModels2 = function(test=FALSE, tausq=.1^2, margVar=.15^2, gamma=-1, 
                             beta0=-1.75, resultType=c("county", "pixel", "EA"), 
                            sampling=c("SRS", "oversamp"), recomputeTruth=TRUE, modelsI=1:15, 
                            produceFigures=FALSE, big=FALSE, printIEvery=50, 
                            maxDataSets=NULL, nsim=10, saveResults=TRUE, loadResults=FALSE, 
                            xtable.args=list(digits=c(0, 1, 1, 1, 1, 0, 1), display=rep("f", 7), auto=TRUE), 
                            tableFormat=c("2", "1"), colScale=c(10^4, 10^5, 100^2, 10^3, 100, 100), 
                            colUnits=c(" ($\\times 10^{-4}$)", " ($\\times 10^{-5}$)", " ($\\times 10^{-4}$)", 
                                       " ($\\times 10^{-3}$)", " ($\\times 10^{-2}$)", " ($\\times 10^{-2}$)"), 
                            colDigits=c(1, 1, 1, 1, 0, 1)) {
  
  # match the arguments with their correct values
  resultType = match.arg(resultType)
  sampling = match.arg(sampling)
  tableFormat = match.arg(tableFormat)
  
  # test to make sure only the naive and direct models are included if big is set to TRUE
  if(!identical(modelsI, 1:2) && big)
    stop("if big is set to TRUE, only naive and direct results can be included")
  
  # load the true superpopulation
  testText = ifelse(test, "Test", "")
  bigText = ifelse(big, "Big", "")
  if(!test)
    load(paste0("simDataMultiBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                "HHoldVar0urbanOverSamplefrac0", bigText, ".RData"))
  else
    load(paste0("simDataMultiBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                "HHoldVar0urbanOverSamplefrac0Test", bigText, ".RData"))
  eaDat = SRSDat$eaDat
  
  if(sampling == "SRS") {
    clustDat = SRSDat
  } else {
    clustDat = overSampDat
  }
  maxDataSets = ifelse(is.null(maxDataSets), length(clustDat$clustDat), maxDataSets)
  
  # allModels = c("naive", "direct", "mercer", "bym", "bymMod", "bymNoUrb", "bymNoUrbMod", "bymNoClust", "bymNoUrbClust", "spde", "spdeNoUrb")
  # allNames = c("Naive", "Direct ", "Mercer et al.", "BYM (no urban/cluster)", "BYM (no urban)", "BYM (no cluster)", "BYM", "SPDE (no urban)", "SPDE")
  # allNamesBinomial = c("Naive Binom.", "Direct Binom.", "Mercer et al. Binom.", "BYM Binom. (no urb/clust)", "BYM Binom. (no urb)", "BYM Binom. (no clust)", "BYM Binom.", "SPDE Binom. (no urb)", "SPDE Binom.")
  # BYM models are in order of complexity: no urban/cluster, no urban, no cluster, full
  allNames = c("Naive", "Direct", "Mercer", "BYM2 1", "BYM2 2", "BYM2 2'", "BYM2 3", "BYM2 4", "BYM2 4'", 
               "BYM2 Pop 1", "BYM2 Pop 2", "BYM2 Pop 2'", "BYM2 Pop 3", "BYM2 Pop 4", "BYM2 Pop 4'", 
               "SPDE 1", "SPDE 2", "SPDE 3", "SPDE 4")
  allNamesBinomial = c("Naive Binom.", "Direct Binom.", "Mercer Binom.", "BYM2 1 Binom.", "BYM2 2 Binom.", "BYM2 2' Binom.", "BYM2 3 Binom.", "BYM2 4 Binom.", "BYM2 4' Binom.", 
                       "BYM2 Pop 1 Binom.", "BYM2 Pop 2 Binom.", "BYM2 Pop 2' Binom.", "BYM2 Pop 3 Binom.", "BYM2 Pop 4 Binom.", "BYM2 Pop 4' Binom.", 
                       "SPDE 1 Binom.", "SPDE 2 Binom.", "SPDE 3 Binom.", "SPDE 4")
  # allNamesBinomial = c("Naive Binom.", "Direct Binom.", "Mercer Binom.", "BYM 1 Binom.", "BYM 2 Binom.", "BYM 2' Binom.", "BYM 3 Binom.", "BYM 4 Binom.", "BYM 4' Binom.", "SPDE 1 Binom.", "SPDE 2 Binom.")
  models = allModels[modelsI]
  
  # this string carries all the information about the run
  runId = paste0("Tausq", round(tausq, 3), "margVar", round(margVar, 3), "gamma", round(gamma, 3), 
                 testText, bigText, sampling, 
                 "models", do.call("paste0", as.list(modelsI)), "nsim", nsim, "MaxDataSetI", maxDataSets)
  
  # compute the results if necessary
  if(!loadResults) {
    
    ## generate the true values at the county level for the 
    ## given settings if these are new settings
    counties = sort(unique(eaDat$admin1))
    if(recomputeTruth){
      truth = getTruth(resultType, eaDat)
      if(resultType == "county") {
        save(truth, file="truthbycounty.RData")
      } else if(resultType == "pixel") {
        save(truth, file="truthbyPixel.RData")
      } else if(resultType == "EA") {
        save(truth, file="truthbyEA.RData")
      }
    } else {
      # load the truth
      load(paste0("truthby", resultType, ".RData"))
    }
    
    ## Now pick which models to include in the results
    
    # load data
    tauText = ifelse(tausq == 0, "0", "0.01")
    testText = ifelse(test, "Test", "")
    if("naive" %in% models || "direct" %in% models) {
      out = load(paste0("resultsDirectNaiveBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                       "HHoldVar0urbanOverSamplefrac0Test", bigText, ".RData"))
      if(sampling == "SRS") {
        directEst = directEstSRS
        naive = naiveSRS
      }
      else {
        directEst = directEstoverSamp
        naive = naiveoverSamp
      }
    }
    if("mercer" %in% models) {
      if(!test)
        load(paste0("resultsMercerBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                                                    "HHoldVar0urbanOverSamplefrac0.RData"))
      else
        load(paste0("resultsMercerBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                                                    "HHoldVar0urbanOverSamplefrac0Test.RData"))
      if(sampling == "SRS")
        mercer = mercerSRS
      else
        mercer = merceroverSamp
    }
    if("bymNoUrb" %in% models) {
      includeUrbanRural = FALSE
      includeCluster = TRUE
      aggregateByPopulation = FALSE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                         includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 100, testText, '.RData'))
      # out = load(paste0("kenyaSpatialDesignResultNewTausq", tauText, "UrbRurFALSEClusterTRUE", testText, ".RData"))
      if(sampling == "SRS")
        designRes$overSampDat = NULL
      else
        designRes$SRSdat = NULL
      designResNoUrb = designRes
    }
    if("bymNoUrbMod" %in% models) {
      out = load(paste0("kenyaSpatialDesignResultNewTausq", tauText, "UrbRurFALSEClusterTRUEdebiased", testText, ".RData"))
      if(sampling == "SRS")
        designRes$overSampDat = NULL
      else
        designRes$SRSdat = NULL
      designResNoUrbMod = designRes
    }
    if("bymNoUrbClust" %in% models) {
      out = load(paste0("kenyaSpatialDesignResultNewTausq", tauText, "UrbRurFALSEClusterFALSE", testText, ".RData"))
      if(sampling == "SRS")
        designRes$overSampDat = NULL
      else
        designRes$SRSdat = NULL
      designResNoUrbClust = designRes
    }
    if("bymNoClust" %in% models) {
      out = load(paste0("kenyaSpatialDesignResultNewTausq", tauText, "UrbRurTRUEClusterFALSE", testText, ".RData"))
      if(sampling == "SRS")
        designRes$overSampDat = NULL
      else
        designRes$SRSdat = NULL
      designResNoClust = designRes
    }
    if("bymMod" %in% models) {
      out = load(paste0("kenyaSpatialDesignResultNewTausq", tauText, "UrbRurTRUEClusterTRUEdebiased", testText, ".RData"))
      if(sampling == "SRS")
        designRes$overSampDat = NULL
      else
        designRes$SRSdat = NULL
      designResMod = designRes
    }
    if("bym" %in% models) {
      out = load(paste0("kenyaSpatialDesignResultNewTausq", tauText, "UrbRurTRUEClusterTRUE", testText, ".RData"))
      if(sampling == "SRS")
        designRes$overSampDat = NULL
      else
        designRes$SRSdat = NULL
    }
    if("spdeNoUrb" %in% models) {
      out = load(paste0("resultsSPDETausq", tauText, "urbanEffectFALSE", testText, ".RData"))
      if(sampling == "SRS")
        spdeNoUrb = spdeSRS
      else
        spdeNoUrb = spdeOverSamp
    }
    if("spde" %in% models) {
      out = load(paste0("resultsSPDETausq", tauText, "urbanEffectTRUE", testText, ".RData"))
      if(sampling == "SRS")
        spde = spdeSRS
      else
        spde = spdeOverSamp
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
    
    # calculate the binomial n (number of children) for each prediction unit
    numChildren = truth$numChildren
    
    # compute scores
    scoresDirect = scoresNaive = scoresMercer = scoresBYMNoUrb = scoresBYM = scoresBYMNoUrbMod = scoresBYMMod = scoresBYMNoUrbClust = scoresBYMNoClust = scoresSPDENoUrb = scoresSPDE = data.frame()
    
    # convert results to the desired aggregation level
    # not including urban effect
    if("bymNoUrb" %in% models) {
      designResNoUrb[[1]]$Q10 = getSubLevelResults(designResNoUrb[[1]]$Q10)
      designResNoUrb[[1]]$Q50 = getSubLevelResults(designResNoUrb[[1]]$Q50)
      designResNoUrb[[1]]$Q90 = getSubLevelResults(designResNoUrb[[1]]$Q90)
      designResNoUrb[[1]]$mean = getSubLevelResults(designResNoUrb[[1]]$mean)
      designResNoUrb[[1]]$stddev = getSubLevelResults(designResNoUrb[[1]]$stddev)
    }
    
    # including urban effect
    if("bym" %in% models) {
      designRes[[1]]$Q10 = getSubLevelResults(designRes[[1]]$Q10)
      designRes[[1]]$Q50 = getSubLevelResults(designRes[[1]]$Q50)
      designRes[[1]]$Q90 = getSubLevelResults(designRes[[1]]$Q90)
      designRes[[1]]$mean = getSubLevelResults(designRes[[1]]$mean)
      designRes[[1]]$stddev = getSubLevelResults(designRes[[1]]$stddev)
    }
    
    # not including urban effect
    if("bymNoUrbClust" %in% models) {
      designResNoUrbClust[[1]]$Q10 = getSubLevelResults(designResNoUrbClust[[1]]$Q10)
      designResNoUrbClust[[1]]$Q50 = getSubLevelResults(designResNoUrbClust[[1]]$Q50)
      designResNoUrbClust[[1]]$Q90 = getSubLevelResults(designResNoUrbClust[[1]]$Q90)
      designResNoUrbClust[[1]]$mean = getSubLevelResults(designResNoUrbClust[[1]]$mean)
      designResNoUrbClust[[1]]$stddev = getSubLevelResults(designResNoUrbClust[[1]]$stddev)
    }
    
    # including urban effect
    if("bymNoClust" %in% models) {
      designResNoClust[[1]]$Q10 = getSubLevelResults(designResNoClust[[1]]$Q10)
      designResNoClust[[1]]$Q50 = getSubLevelResults(designResNoClust[[1]]$Q50)
      designResNoClust[[1]]$Q90 = getSubLevelResults(designResNoClust[[1]]$Q90)
      designResNoClust[[1]]$mean = getSubLevelResults(designResNoClust[[1]]$mean)
      designResNoClust[[1]]$stddev = getSubLevelResults(designResNoClust[[1]]$stddev)
    }
    
    # not including urban effect, modified to be debiased using marginal rather than conditional effect as prediction
    if("bymNoUrbMod" %in% models) {
      designResNoUrbMod[[1]]$Q10 = getSubLevelResults(designResNoUrbMod[[1]]$Q10)
      designResNoUrbMod[[1]]$Q50 = getSubLevelResults(designResNoUrbMod[[1]]$Q50)
      designResNoUrbMod[[1]]$Q90 = getSubLevelResults(designResNoUrbMod[[1]]$Q90)
      designResNoUrbMod[[1]]$mean = getSubLevelResults(designResNoUrbMod[[1]]$mean)
      designResNoUrbMod[[1]]$stddev = getSubLevelResults(designResNoUrbMod[[1]]$stddev)
    }
    
    # including urban effect, modified to be debiased using marginal rather than conditional effect as prediction
    if("bymMod" %in% models) {
      designResMod[[1]]$Q10 = getSubLevelResults(designResMod[[1]]$Q10)
      designResMod[[1]]$Q50 = getSubLevelResults(designResMod[[1]]$Q50)
      designResMod[[1]]$Q90 = getSubLevelResults(designResMod[[1]]$Q90)
      designResMod[[1]]$mean = getSubLevelResults(designResMod[[1]]$mean)
      designResMod[[1]]$stddev = getSubLevelResults(designResMod[[1]]$stddev)
    }
    
    for(i in c(1:maxDataSets)) { # for problem fitting mercerSRS for SRS sampling, tausq=0
      # for(i in 1:100) {
      if((i %% printIEvery == 0) || (i == 1))
        print(i)
      resultName = paste0(resultType, "Results")
      if(resultType == "EA")
        resultName = "eaResults"
      
      # convert results to the desired aggregation level
      if("direct" %in% models)
        directEsti = getSubLevelResults(directEst[[i]])
      if("naive" %in% models)
        naivei = getSubLevelResults(naive[[i]])
      if("mercer" %in% models)
        merceri = getSubLevelResults(mercer[[i]])
      if(resultType != "county") {
        if("spdeNoUrb" %in% models)
          spdeNoUrbi = spdeNoUrb[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
        if("spde" %in% models)
          spdei = spde[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
        if("direct" %in% models)
          directEsti = getSubLevelResults(directEst[[i]])
      } else {
        if("spdeNoUrb" %in% models)
          spdeNoUrbi = spdeNoUrb[[resultName]][[i]]
        if("spde" %in% models)
          spdei = spde[[resultName]][[i]]
      }
      
      if(resultType == "EA") {
        # set first row of spde results to be the EA index
        if("spdeNoUrb" %in% models) {
          spdeNoUrbi[[resultType]] = 1:nrow(spdeNoUrbi)
          
          whichName = which(names(spdeNoUrbi) == "EA")
          spdeNoUrbi = cbind(spdeNoUrbi[,whichName], spdeNoUrbi[,-whichName])
          
          names(spdeNoUrbi)[1] = "EA"
        }
        if("spde" %in% models) {
          spdei[[resultType]] = 1:nrow(spdei)
          
          whichName = which(names(spdei) == "EA")
          spdei = cbind(spdei[,whichName], spdei[,-whichName])
          
          names(spdei)[1] = "EA"
        }
      }
      
      # for spde results, modify the name of the results
      # modify result row and column table names according to aggregation level
      if(resultType == "county") {
        # without urban effect:
        if("spdeNoUrb" %in% models) {
          whichName = which(names(spdeNoUrbi) == "admin1")
          names(spdeNoUrbi)[whichName] = resultType
        }
        
        if("spde" %in% models) {
          # with urban effect:
          whichName = which(names(spdei) == "admin1")
          names(spdei)[whichName] = resultType
        }
      }
      
      if(resultType == "pixel") {
        # set first row of spde results to be the pixel index
        if("spdeNoUrb" %in% models) {
          spdeNoUrbi[[resultType]] = truth$pixel
          
          whichName = which(names(spdeNoUrbi) == "pixel")
          spdeNoUrbi = cbind(spdeNoUrbi[,whichName], spdeNoUrbi[,-whichName])
          
          names(spdeNoUrbi)[1] = "pixel"
        }
        if("spde" %in% models) {
          spdei[[resultType]] = truth$pixel
          
          whichName = which(names(spdei) == "pixel")
          spdei = cbind(spdei[,whichName], spdei[,-whichName])
          
          names(spdei)[1] = "pixel"
        }
      }
      
      # change names of table variables in spde model with no urban effect to reflect that
      if("spdeNoUrb" %in% models)
        names(spdeNoUrbi)[2:6] = paste0(names(spdeNoUrbi)[2:6], "NoUrb")
      
      if("direct" %in% models) {
        allres = merge(truth, directEsti, by=resultType)
        colnames(allres) = c(resultType, "truth", paste(colnames(allres)[3:8], "direct", sep=""))
      } else {
        stop("direct estimates must be included at this point in order to name the estimate table columns")
      }
      if("naive" %in% models)
        allres = merge(allres, naivei, by=resultType)
      if("mercer" %in% models)
        allres = merge(allres, merceri, by=resultType)
      if("spdeNoUrb" %in% models)
        allres = merge(allres, spdeNoUrbi, by=resultType)
      if("spde" %in% models)
        allres = merge(allres, spdei, by=resultType)
      
      # set whether or not to calculate scores on logit scale depending on result type
      thisTruth = allres$truth
      
      # if(is.null(useLogit)) {
      #   useLogit = FALSE
      #   if(resultType != "EA" && resultType != "pixel") {
      #     thisTruth = logit(thisTruth)
      #     useLogit=TRUE
      #   }
      # } else if(useLogit == TRUE) {
      #   thisTruth = logit(thisTruth)
      # }
      
      if("direct" %in% models) {
        my.scoresdirect = getScores(thisTruth, numChildren, allres$logit.estdirect, allres$var.estdirect, nsim=nsim)
        scoresDirect <- rbind(scoresDirect,
                              cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresdirect))
      }
      if("naive" %in% models) {
        my.scoresnaive = getScores(thisTruth, numChildren, allres$logit.est, allres$var.est, nsim=nsim)
        scoresNaive <- rbind(scoresNaive,
                             cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresnaive))
      }
      if("mercer" %in% models) {
        my.scoresmercer = getScores(thisTruth, numChildren, allres$logit.est.mercer, allres$var.est.mercer, nsim=nsim)
        scoresMercer <- rbind(scoresMercer,
                              cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresmercer))
      }
      if("bymNoUrbClust" %in% models) {
        my.scoresbymNoUrbClust = getScores(thisTruth, numChildren, designResNoUrbClust[[1]]$mean[,i], (designResNoUrbClust[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoUrbClust <- rbind(scoresBYMNoUrbClust,
                                     cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoUrbClust))
      }
      if("bymNoUrb" %in% models) {
        my.scoresbymNoUrb = getScores(thisTruth, numChildren, designResNoUrb[[1]]$mean[,i], (designResNoUrb[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoUrb <- rbind(scoresBYMNoUrb,
                                cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoUrb))
      }
      if("bymNoUrbMod" %in% models) {
        my.scoresbymNoUrbMod = getScores(thisTruth, numChildren, designResNoUrbMod[[1]]$mean[,i], (designResNoUrbMod[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoUrbMod <- rbind(scoresBYMNoUrbMod,
                                cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoUrbMod))
      }
      if("bymNoClust" %in% models) {
        my.scoresbymNoClust = getScores(thisTruth, numChildren, designResNoClust[[1]]$mean[,i], (designResNoClust[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoClust <- rbind(scoresBYMNoClust,
                                  cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoClust))
      }
      if("bym" %in% models) {
        my.scoresbym = getScores(thisTruth, numChildren, designRes[[1]]$mean[,i], (designRes[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYM <- rbind(scoresBYM,
                           cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbym))
      }
      if("bymMod" %in% models) {
        my.scoresbymMod = getScores(thisTruth, numChildren, designResMod[[1]]$mean[,i], (designResMod[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMMod <- rbind(scoresBYMMod,
                           cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymMod))
      }
      if("spdeNoUrb" %in% models) {
        stop("determine if the spde code should compute all of these directly")
        my.biasspdeNoUrb = bias(thisTruth, allres$logit.est.spdeNoUrb, logit=useLogit, my.var=allres$var.est.spdeNoUrb)
        my.msespdeNoUrb = mse(thisTruth, allres$logit.est.spdeNoUrb, logit=useLogit, my.var=allres$var.est.spdeNoUrb)
        my.dssspdeNoUrb = dss(thisTruth, allres$logit.est.spdeNoUrb, allres$var.est.spdeNoUrb)
        # the below line used to be commented for some reason
        my.crpsspdeNoUrb = crpsNormal(thisTruth, allres$logit.est.spdeNoUrb, allres$var.est.spdeNoUrb, logit=useLogit, n=numChildren)
        my.crpsspdeNoUrb = allres$crps.spdeNoUrb
        my.coveragespdeNoUrb = coverage(thisTruth, allres$lower.spdeNoUrb, allres$upper.spdeNoUrb, logit=useLogit)
        my.lengthspdeNoUrb = intervalWidth(allres$lower.spdeNoUrb, allres$upper.spdeNoUrb, logit=useLogit)
        scoresSPDENoUrb <- rbind(scoresSPDENoUrb,
                                 cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresspdeNoUrb))
      }
      if("spde" %in% models) {
        stop("determine if the spde code should compute all of these directly")
        my.biasspde = bias(thisTruth, allres$logit.est.spde, logit=useLogit, my.var=allres$var.est.spde)
        my.msespde = mse(thisTruth, allres$logit.est.spde, logit=useLogit, my.var=allres$var.est.spde)
        my.dssspde = dss(thisTruth, allres$logit.est.spde, allres$var.est.spde)
        # the below line needs to be commented for some reason
        my.crpsspde = crpsNormal(thisTruth, allres$logit.est.spde, allres$var.est.spde, logit=useLogit, n=numChildren)
        my.crpsspde = allres$crps.spde
        my.coveragespde = coverage(thisTruth, allres$lower.spde, allres$upper.spde, logit=useLogit)
        my.lengthspde = intervalWidth(allres$lower.spde, allres$upper.spde, logit=useLogit)
        scoresSPDE <- rbind(scoresSPDE,
                            cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresspde))
      }
    }
  }
  else {
    # in this case, we have already computed the results so just load them into the environment
    load(paste0("scores", runId, ".RData"))
    
    # also, don't resave these results that we've already saved
    saveResults = FALSE
  }
  
  # final table: 
  if("naive" %in% models)
    naive = apply(scoresNaive[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("direct" %in% models)
    direct = apply(scoresDirect[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("mercer" %in% models)
    mercer = apply(scoresMercer[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("bymNoUrbClust" %in% models)
    bymNoUrbClust = apply(scoresBYMNoUrbClust[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("bymNoUrb" %in% models)
    bymNoUrb = apply(scoresBYMNoUrb[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("bymNoUrbMod" %in% models)
    bymNoUrbMod = apply(scoresBYMNoUrbMod[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("bymNoClust" %in% models)
    bymNoClust = apply(scoresBYMNoClust[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("bym" %in% models)
    bym = apply(scoresBYM[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("bymMod" %in% models)
    bymMod = apply(scoresBYMMod[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("spdeNoUrb" %in% models)
    spdeNoUrb = apply(scoresSPDENoUrb[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("spde" %in% models)
    spde = apply(scoresSPDE[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  idx = 1:9
  tab = c()
  if("naive" %in% models)
    tab = rbind(tab, c(naive[idx]))
  if("direct" %in% models)
    tab = rbind(tab, c(direct[idx]))
  if("mercer" %in% models)
    tab = rbind(tab, c(mercer[idx]))
  if("bymNoUrbClust" %in% models)
    tab = rbind(tab, c(bymNoUrbClust[idx]))
  if("bymNoUrb" %in% models)
    tab = rbind(tab, c(bymNoUrb[idx]))
  if("bymNoUrbMod" %in% models)
    tab = rbind(tab, c(bymNoUrbMod[idx]))
  if("bymNoClust" %in% models)
    tab = rbind(tab, c(bymNoClust[idx]))
  if("bym" %in% models)
    tab = rbind(tab, c(bym[idx]))
  if("bymMod" %in% models)
    tab = rbind(tab, c(bymMod[idx]))
  if("spdeNoUrb" %in% models)
    tab = rbind(tab, c(spdeNoUrb[idx]))
  if("spde" %in% models)
    tab = rbind(tab, c(spde[idx]))
  colnames(tab) = c("Bias", "Var", "MSE", "CRPS", "CRPS Binom.", "80\\% Cvg", "80\\% Cvg Binom.", "CI Width", "CI Width Binom.")
  finalNames = allNames[modelsI]
  finalNamesBinomial = allNamesBinomial[modelsI]
  
  if(tableFormat == "1")
    rownames(tab) = finalNames
  else {
    # separate out the binomial and non-binomial models into separate rows
    binomialColumns = c(1:3, 5, 7, 9)
    otherColumns = c(1:3, 4, 6, 8)
    binomialTable = tab[, binomialColumns]
    otherTable = tab[, otherColumns]
    tab = otherTable[numeric(0),]
    thisFinalNames = c()
    
    for(i in 1:length(modelsI)) {
      tab = rbind(tab, 
                  otherTable[i,], 
                  binomialTable[i,])
      thisFinalNames = c(thisFinalNames, finalNames[i], finalNamesBinomial[i])
    }
    rownames(tab) = thisFinalNames
  }
  
  # round the columns of tab and modify the column names to include the scale
  for(i in 1:ncol(tab)) {
    tab[,i] = as.numeric(print(round(tab[,i] * colScale[i], digits=colDigits[i])))
    colnames(tab)[i] = paste0(colnames(tab)[i], colUnits[i])
  }
  
  print(do.call("xtable", c(list(tab), xtable.args)), 
        include.colnames=TRUE,
        hline.after=0, 
        math.style.exponents=TRUE, 
        sanitize.text.function=function(x){x})
  
  runId = paste0("Tausq", round(tausq, 3), testText, bigText, sampling, 
                 "models", do.call("paste0", as.list(modelsI)), "nsim", nsim, "MaxDataSetI", maxDataSets)
  if(saveResults) {
    # first collect all the results. Save everything except for the postprocessing arguments: 
    # produceFigures, digits
    objectNames = ls()
    objectNames = objectNames[-match(c("produceFigures", "xtable.args", "tableFormat", "colScale", 
                                       "colUnits", "colDigits"), objectNames)]
    save(list=objectNames, file=paste0("scores", runId, ".RData"))
  }
  
  if(produceFigures) {
    # compare all five
    pdf(paste0("figures/biasbyregion", runId, ".pdf"), width=20, height=12)
    par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
    boxplot(bias~region, data=scoresDirect, at=seq(-1, 275, by=6), 
            col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
    boxplot(bias~region, data=scoresNaive, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
    boxplot(bias~region, data=scoresMercer, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
    boxplot(bias~region, data=scoresBYM, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
    boxplot(bias~region, data=scoresSPDE, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
    # axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirect$region[1:47])
    # axis(1, at=seq(0.5, 276.5, by=6), labels=scoresDirect$region[1:47])
    legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
           fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
    abline(h=0, lwd=2, col=2)
    dev.off()
    
    pdf(paste0("figures/crpsbyregion", runId, ".pdf"), width=20, height=12)
    par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
    boxplot(crps~region, data=scoresDirect, at=seq(-1, 275, by=6), 
            col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
    if("naive" %in% models)
      boxplot(crps~region, data=scoresNaive, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
    if("mercer" %in% models)
      boxplot(crps~region, data=scoresMercer, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
    if("BYM" %in% models)
      boxplot(crps~region, data=scoresBYM, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
    if("spde" %in% models)
      boxplot(crps~region, data=scoresSPDE, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
    # axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirect$region[1:47])
    legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
           fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
    dev.off()
  }
  
  tab
}

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
# logitScale: whether or not to produce scoring rule results at the logit or count scale. 
#           By default this is to use loaded scale results unless the resultType is EA or pixel
# printIEvery: how often to print progress
runCompareModels = function(test=FALSE, tausq=.1^2, resultType=c("county", "pixel", "EA"), 
                            sampling=c("SRS", "oversamp"), recomputeTruth=TRUE, modelsI=1:3, 
                            produceFigures=FALSE, logitScale=NULL, big=FALSE, printIEvery=50) {
  
  # match the arguments with their correct values
  resultType = match.arg(resultType)
  sampling = match.arg(sampling)
  useLogit = logitScale
  
  # test to make sure only the naive and direct models are included if big is set to TRUE
  if(!identical(modelsI, 1:2) && big)
    stop("if big is set to TRUE, only naive and direct results can be included")
  
  testText = ifelse(test, "Test", "")
  bigText = ifelse(big, "Big", "")
  if(tausq == .1^2 || tausq == .01) {
    out = load(paste0("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOverSamplefrac0", testText, bigText, ".RData"))
  }
  else {
    if(tausq != 0)
      stop("tausq can only be equal to .1^2 or 0")
    
    out = load(paste0("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOverSamplefrac0", testText, bigText, ".RData"))
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
      childrenPerCounty = truthbycounty
      
      for(i in 1:47){
        super = eaDat[eaDat$admin1 == regions[i],]
        childrenPerCounty[i] = sum(super$numChildren)
        truthbycounty[i] <- sum(super$died)/childrenPerCounty[i]
      }
      truth = data.frame(admin1=regions, truth=truthbycounty, numChildren=childrenPerCounty)
      save(truth, file="truthbycounty.RData")
    } else if(resultType == "pixel") {
      counties = sort(unique(eaDat$admin1))
      eaToPixel = eaDat$pixelI
      childrenPerPixel = tapply(eaDat$numChildren, list(pixel=eaDat$pixelI), sum)
      urbanPixel = tapply(eaDat$urban, list(pixel=eaDat$pixelI), function(x) {mean(x[1])})
      deathsPerPixel = tapply(eaDat$died, list(pixel=eaDat$pixelI), sum)
      regions = names(childrenPerPixel) # these are the pixels with enumeration areas in them
      pixelToAdmin = match(popGrid$admin1[as.numeric(regions)], counties)
      
      truth = data.frame(pixel=regions, truth=deathsPerPixel / childrenPerPixel, countyI=pixelToAdmin, urban=urbanPixel, numChildren=childrenPerPixel)
      save(truth, file="truthbyPixel.RData")
    } else if(resultType == "EA") {
      truth = data.frame(EA = 1:nrow(eaDat), truth = eaDat$died/25, urban=eaDat$urban, numChildren=eaDat$numChildren)
      save(truth, file="truthbyEA.RData")
    }
  } else {
    # load the truth
    load(paste0("truthby", resultType, ".RData"))
  }
  
  ## Now pick which models to include in the results
  allModels = c("naive", "direct", "mercer", "bym", "bymNoUrb", "bymNoClust", "bymNoUrbClust", "spde", "spdeNoUrb")
  allNames = c("Naive", "Direct estimates", "Mercer et al.", "Model-based BYM", "Model-based BYM (no urban effect)", "Model-based BYM (no cluster effect)", "Model-based BYM (no urban or cluster effect)", "SPDE", "SPDE (no urban effect)")
  models = allModels[modelsI]
  
  # load data
  if(tausq == .1^2 || tausq == 0.01) {
    if(test) {
      if("naive" %in% models || "direct" %in% models)
        out = load(paste0("resultsDirectNaiveTausq0.01Test", bigText, ".RData")) # this is the only case that uses the big dataset
      if("mercer" %in% models)
        out = load("resultsMercerTausq0.01test.RData")
      if("bymNoUrb" %in% models) {
        out = load("kenyaSpatialDesignResultNewTausq0.01UrbRurFALSEClusterTRUETest.RData")
        designResNoUrb = designRes
      }
      if("bym" %in% models)
        out = load("kenyaSpatialDesignResultNewTausq0.01UrbRurTRUEClusterTRUETest.RData")
      if("bymNoUrbClust" %in% models) {
        out = load("kenyaSpatialDesignResultNewTausq0.01UrbRurFALSEClusterFALSETest.RData")
        designResNoUrb = designRes
      }
      if("bymNoClust" %in% models)
        out = load("kenyaSpatialDesignResultNewTausq0.01UrbRurTRUEClusterFALSETest.RData")
      if("spdeNoUrb" %in% models) {
        out = load("resultsSPDETausq0.01urbanEffectFALSETest.RData")
        spdeSRSNoUrb = spdeSRS
        spdeOverSampNoUrb = spdeOverSamp
      }
      if("spde" %in% models)
        out = load("resultsSPDETausq0.01urbanEffectTRUETest.RData")
    } else {
      if("naive" %in% models || "direct" %in% models)
        out = load(paste0("resultsDirectNaiveTausq0.01", bigText, ".RData"))
      if("mercer" %in% models)
        out = load("resultsMercerTausq0.01.RData")
      if("bymNoUrb" %in% models) {
        out = load("kenyaSpatialDesignResultNewTausq0.01UrbRurFALSEClusterTRUE.RData")
        designResNoUrb = designRes
      }
      if("bym" %in% models)
        out = load("kenyaSpatialDesignResultNewTausq0.01UrbRurTRUEClusterTRUE.RData")
      if("bymNoUrbClust" %in% models) {
        out = load("kenyaSpatialDesignResultNewTausq0.01UrbRurFALSEClusterFALSE.RData")
        designResNoUrb = designRes
      }
      if("bymNoClust" %in% models)
        out = load("kenyaSpatialDesignResultNewTausq0.01UrbRurTRUEClusterFALSE.RData")
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
        out = load(paste0("resultsDirectNaiveTausq0Test", bigText, ".RData")) # this is the only case that uses the big dataset
      if("mercer" %in% models)
        out = load("resultsMercerTausq0test.RData")
      if("bymNoUrb" %in% models) {
        out = load("kenyaSpatialDesignResultNewTausq0UrbRurFATestLSE.RData")
        designResNoUrb = designRes
      }
      if("bym" %in% models)
        out = load("kenyaSpatialDesignResultNewTausq0UrbRurTRUETest.RData")
      if("spdeNoUrb" %in% models) {
        out = load("resultsSPDETausq0urbanEffectFALSETest.RData")
        spdeSRSNoUrb = spdeSRS
        spdeOverSampNoUrb = spdeOverSamp
      }
      if("spde" %in% models)
        out = load("resultsSPDETausq0urbanEffectTRUETest.RData")
    } else {
      if("naive" %in% models || "direct" %in% models)
        out = load(paste0("resultsDirectNaiveTausq0", bigText, ".RData")) # this is the only case that uses the big dataset
      if("mercer" %in% models)
        out = load("resultsMercerTausq0.RData")
      if("bymNoUrb" %in% models) {
        out = load("kenyaSpatialDesignResultNewTausq0UrbRurFALSE.RData")
        designResNoUrb = designRes
      }
      if("bym" %in% models)
        out = load("kenyaSpatialDesignResultNewTausq0UrbRurTRUE.RData")
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
  
  # calculate the binomial n (number of children) for each prediction unit
  numChildren = truth$numChildren
  
  # compute scores
  scoresDirectSRS = scoresNaiveSRS = scoresMercerSRS = scoresBYMNoUrbSRS = scoresBYMSRS = scoresSPDENoUrbSRS = scoresSPDESRS = data.frame()
  scoresDirectoverSamp = scoresNaiveoverSamp = scoresMerceroverSamp = scoresBYMNoUrboverSamp = scoresBYMoverSamp = scoresSPDENoUrboverSamp = scoresSPDEoverSamp = data.frame()
  # scoresDirectSRS = scoresNaiveSRS = scoresMercerSRS = scoresBYMSRS = data.frame()
  # scoresDirectoverSamp = scoresNaiveoverSamp = scoresMerceroverSamp = scoresBYMoverSamp = data.frame()
  
  # convert results to the desired aggregation level
  # not including urban effect
  if("bymNoUrb" %in% models) {
    designResNoUrb[[1]]$Q10 = getSubLevelResults(designResNoUrb[[1]]$Q10)
    designResNoUrb[[1]]$Q50 = getSubLevelResults(designResNoUrb[[1]]$Q50)
    designResNoUrb[[1]]$Q90 = getSubLevelResults(designResNoUrb[[1]]$Q90)
    designResNoUrb[[1]]$mean = getSubLevelResults(designResNoUrb[[1]]$mean)
    designResNoUrb[[1]]$stddev = getSubLevelResults(designResNoUrb[[1]]$stddev)
    designResNoUrb$overSampDat$Q10 = getSubLevelResults(designResNoUrb$overSampDat$Q10)
    designResNoUrb$overSampDat$Q50 = getSubLevelResults(designResNoUrb$overSampDat$Q50)
    designResNoUrb$overSampDat$Q90 = getSubLevelResults(designResNoUrb$overSampDat$Q90)
    designResNoUrb$overSampDat$mean = getSubLevelResults(designResNoUrb$overSampDat$mean)
    designResNoUrb$overSampDat$stddev = getSubLevelResults(designResNoUrb$overSampDat$stddev)
  }
  
  # including urban effect
  if("bym" %in% models) {
    designRes[[1]]$Q10 = getSubLevelResults(designRes[[1]]$Q10)
    designRes[[1]]$Q50 = getSubLevelResults(designRes[[1]]$Q50)
    designRes[[1]]$Q90 = getSubLevelResults(designRes[[1]]$Q90)
    designRes[[1]]$mean = getSubLevelResults(designRes[[1]]$mean)
    designRes[[1]]$stddev = getSubLevelResults(designRes[[1]]$stddev)
    designRes$overSampDat$Q10 = getSubLevelResults(designRes$overSampDat$Q10)
    designRes$overSampDat$Q50 = getSubLevelResults(designRes$overSampDat$Q50)
    designRes$overSampDat$Q90 = getSubLevelResults(designRes$overSampDat$Q90)
    designRes$overSampDat$mean = getSubLevelResults(designRes$overSampDat$mean)
    designRes$overSampDat$stddev = getSubLevelResults(designRes$overSampDat$stddev)
  }
  
  for(i in c(1:length(clustDat$clustDat))) { # for problem fitting mercerSRS for SRS sampling, tausq=0
    # for(i in 1:100) {
    if((i %% printIEvery == 0) || (i == 1))
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
      if("naive" %in% models)
        allresoverSamp = merge(allresoverSamp, naiveoverSampi, by=resultType)
      if("mercer" %in% models)
        allresoverSamp = merge(allresoverSamp, merceroverSampi, by=resultType)
      if("spdeNoUrb" %in% models)
        allresoverSamp = merge(allresoverSamp, spdeOverSampNoUrbi, by=resultType)
      if("spde" %in% models)
        allresoverSamp = merge(allresoverSamp, spdeOverSampi, by=resultType)
    }
    
    # set whether or not to calculate scores on logit scale depending on result type
    if(sampling == "SRS")
      thisTruth = allresSRS$truth
    else
      thisTruth = allresoverSamp$truth
    
    if(is.null(useLogit)) {
      useLogit = FALSE
      if(resultType != "EA" && resultType != "pixel") {
        thisTruth = logit(thisTruth)
        useLogit=TRUE
      }
    } else if(useLogit == TRUE) {
      thisTruth = logit(thisTruth)
    }
    
    # SRS setting 
    if(sampling == "SRS") {
      if("direct" %in% models) {
        my.biasSRSdirect = bias(thisTruth, allresSRS$logit.estdirect, logit=useLogit, my.var=allresSRS$var.estdirect)
        my.mseSRSdirect = mse(thisTruth, allresSRS$logit.estdirect, logit=useLogit, my.var=allresSRS$var.estdirect)
        my.dssSRSdirect = dss(thisTruth, allresSRS$logit.estdirect, allresSRS$var.estdirect)
        my.crpsSRSdirect = crpsNormal(thisTruth, allresSRS$logit.estdirect, allresSRS$var.estdirect, logit=useLogit, n=numChildren)
        my.coverageSRSdirect = coverage(thisTruth, allresSRS$lowerdirect, allresSRS$upperdirect, doLogit=useLogit)
        my.lengthSRSdirect = intervalWidth(allresSRS$lowerdirect, allresSRS$upperdirect, logit=useLogit)
      }
      
      if("naive" %in% models) {
        my.biasSRSnaive = bias(thisTruth, allresSRS$logit.est, logit=useLogit, my.var=allresSRS$var.est)
        my.mseSRSnaive = mse(thisTruth, allresSRS$logit.est, logit=useLogit, my.var=allresSRS$var.est)
        my.dssSRSnaive = dss(thisTruth, allresSRS$logit.est, allresSRS$var.est)
        my.crpsSRSnaive = crpsNormal(thisTruth, allresSRS$logit.est, allresSRS$var.est, logit=useLogit, n=numChildren)
        my.coverageSRSnaive = coverage(thisTruth, allresSRS$lower, allresSRS$upper, doLogit=useLogit)
        my.lengthSRSnaive = intervalWidth(allresSRS$lower, allresSRS$upper, logit=useLogit)
      }
      
      if("mercer" %in% models) {
        my.biasSRSmercer = bias(thisTruth, allresSRS$logit.est.mercer, logit=useLogit, my.var=allresSRS$var.est.mercer)
        my.mseSRSmercer = mse(thisTruth, allresSRS$logit.est.mercer, logit=useLogit, my.var=allresSRS$var.est.mercer)
        my.dssSRSmercer = dss(thisTruth, allresSRS$logit.est.mercer, allresSRS$var.est.mercer)
        my.crpsSRSmercer = crpsNormal(thisTruth, allresSRS$logit.est.mercer, allresSRS$var.est.mercer, logit=useLogit, n=numChildren)
        my.coverageSRSmercer = coverage(thisTruth, allresSRS$lower.mercer, allresSRS$upper.mercer, doLogit=useLogit)
        my.lengthSRSmercer = intervalWidth(allresSRS$lower.mercer, allresSRS$upper.mercer, logit=useLogit)
      }
      
      if("bymNoUrb" %in% models) {
        my.biasSRSbymNoUrb = bias(thisTruth, designResNoUrb[[1]]$mean[,i], logit=useLogit, my.var=(designResNoUrb[[1]]$stddev[,i])^2)
        my.mseSRSbymNoUrb = mse(thisTruth, designResNoUrb[[1]]$mean[,i], logit=useLogit, my.var=(designResNoUrb[[1]]$stddev[,i])^2)
        my.dssSRSbymNoUrb = dss(thisTruth, designResNoUrb[[1]]$mean[,i], (designResNoUrb[[1]]$stddev[,i])^2)
        my.crpsSRSbymNoUrb = crpsNormal(thisTruth, designResNoUrb[[1]]$mean[,i], (designResNoUrb[[1]]$stddev[,i])^2, logit=useLogit, n=numChildren)
        my.coverageSRSbymNoUrb = coverage(thisTruth, designResNoUrb[[1]]$Q10[,i],designResNoUrb[[1]]$Q90[,i], doLogit=useLogit)
        my.lengthSRSbymNoUrb = intervalWidth(designResNoUrb[[1]]$Q10[,i], designResNoUrb[[1]]$Q90[,i], logit=useLogit)
      }
      
      if("bym" %in% models) {
        my.biasSRSbym = bias(thisTruth, designRes[[1]]$mean[,i], logit=useLogit, my.var=(designRes[[1]]$stddev[,i])^2)
        my.mseSRSbym = mse(thisTruth, designRes[[1]]$mean[,i], logit=useLogit, my.var=(designRes[[1]]$stddev[,i])^2)
        my.dssSRSbym = dss(thisTruth, designRes[[1]]$mean[,i], (designRes[[1]]$stddev[,i])^2)
        my.crpsSRSbym = crpsNormal(thisTruth, designRes[[1]]$mean[,i], (designRes[[1]]$stddev[,i])^2, logit=useLogit, n=numChildren)
        my.coverageSRSbym = coverage(thisTruth, designRes[[1]]$Q10[,i],designRes[[1]]$Q90[,i], doLogit=useLogit)
        my.lengthSRSbym = intervalWidth(designRes[[1]]$Q10[,i], designRes[[1]]$Q90[,i], logit=useLogit)
      }
      
      if("spdeNoUrb" %in% models) {
        my.biasSRSspdeNoUrb = bias(thisTruth, allresSRS$logit.est.spdeNoUrb, logit=useLogit, my.var=allresSRS$var.est.spdeNoUrb)
        my.mseSRSspdeNoUrb = mse(thisTruth, allresSRS$logit.est.spdeNoUrb, logit=useLogit, my.var=allresSRS$var.est.spdeNoUrb)
        my.dssSRSspdeNoUrb = dss(thisTruth, allresSRS$logit.est.spdeNoUrb, allresSRS$var.est.spdeNoUrb)
        # the below line used to be commented for some reason
        my.crpsSRSspdeNoUrb = crpsNormal(thisTruth, allresSRS$logit.est.spdeNoUrb, allresSRS$var.est.spdeNoUrb, logit=useLogit, n=numChildren)
        my.crpsSRSspdeNoUrb = allresSRS$crps.spdeNoUrb
        my.coverageSRSspdeNoUrb = coverage(thisTruth, allresSRS$lower.spdeNoUrb, allresSRS$upper.spdeNoUrb, doLogit=useLogit)
        my.lengthSRSspdeNoUrb = intervalWidth(allresSRS$lower.spdeNoUrb, allresSRS$upper.spdeNoUrb, logit=useLogit)
      }
      
      if("spde" %in% models) {
        my.biasSRSspde = bias(thisTruth, allresSRS$logit.est.spde, logit=useLogit, my.var=allresSRS$var.est.spde)
        my.mseSRSspde = mse(thisTruth, allresSRS$logit.est.spde, logit=useLogit, my.var=allresSRS$var.est.spde)
        my.dssSRSspde = dss(thisTruth, allresSRS$logit.est.spde, allresSRS$var.est.spde)
        # the below line needs to be commented for some reason
        my.crpsSRSspde = crpsNormal(thisTruth, allresSRS$logit.est.spde, allresSRS$var.est.spde, logit=useLogit, n=numChildren)
        my.crpsSRSspde = allresSRS$crps.spde
        my.coverageSRSspde = coverage(thisTruth, allresSRS$lower.spde, allresSRS$upper.spde, doLogit=useLogit)
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
                                              var=mean((designResNoUrb[[1]]$stddev[,i])^2),
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
                                         var=mean((designRes[[1]]$stddev[,i])^2),
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
        my.crpsoverSampdirect = crpsNormal(thisTruth, allresoverSamp$logit.estdirect, allresoverSamp$var.estdirect, logit=useLogit, n=numChildren)
        my.coverageoverSampdirect = coverage(thisTruth, allresoverSamp$upperdirect, allresoverSamp$lowerdirect, doLogit=useLogit)
        my.lengthoverSampdirect = intervalWidth(allresoverSamp$lowerdirect, allresoverSamp$upperdirect, logit=useLogit)
      }
      
      if("naive" %in% models) {
        my.biasoverSampnaive = bias(thisTruth, allresoverSamp$logit.est, logit=useLogit, my.var=allresoverSamp$var.est)
        my.mseoverSampnaive = mse(thisTruth, allresoverSamp$logit.est, logit=useLogit, my.var=allresoverSamp$var.est)
        my.dssoverSampnaive = dss(thisTruth, allresoverSamp$logit.est, allresoverSamp$var.est)
        my.crpsoverSampnaive = crpsNormal(thisTruth, allresoverSamp$logit.est, allresoverSamp$var.est, logit=useLogit, n=numChildren)
        my.coverageoverSampnaive = coverage(thisTruth, allresoverSamp$upper, allresoverSamp$lower, doLogit=useLogit)
        my.lengthoverSampnaive = intervalWidth(allresoverSamp$lower, allresoverSamp$upper, logit=useLogit)
      }
      
      if("mercer" %in% models) {
        my.biasoverSampmercer = bias(thisTruth, allresoverSamp$logit.est.mercer, logit=useLogit, my.var=allresoverSamp$var.est.mercer)
        my.mseoverSampmercer = mse(thisTruth, allresoverSamp$logit.est.mercer, logit=useLogit, my.var=allresoverSamp$var.est.mercer)
        my.dssoverSampmercer = dss(thisTruth, allresoverSamp$logit.est.mercer, allresoverSamp$var.est.mercer)
        my.crpsoverSampmercer = crpsNormal(thisTruth, allresoverSamp$logit.est.mercer, allresoverSamp$var.est.mercer, logit=useLogit, n=numChildren)
        my.coverageoverSampmercer = coverage(thisTruth, allresoverSamp$lower.mercer, allresoverSamp$upper.mercer, doLogit=useLogit)
        my.lengthoverSampmercer = intervalWidth(allresoverSamp$lower.mercer, allresoverSamp$upper.mercer, logit=useLogit)
      }
      
      if("bymNoUrb" %in% models) {
        my.biasoverSampbymNoUrb = bias(thisTruth, designResNoUrb$overSampDat$mean[,i], logit=useLogit, my.var=(designResNoUrb$overSampDat$stddev[,i])^2)
        my.mseoverSampbymNoUrb = mse(thisTruth, designResNoUrb$overSampDat$mean[,i], logit=useLogit, my.var=(designResNoUrb$overSampDat$stddev[,i])^2)
        my.dssoverSampbymNoUrb = dss(thisTruth, designResNoUrb$overSampDat$mean[,i], (designResNoUrb$overSampDat$stddev[,i])^2)
        my.crpsoverSampbymNoUrb = crpsNormal(thisTruth, designResNoUrb$overSampDat$mean[,i], (designResNoUrb$overSampDat$stddev[,i])^2, logit=useLogit, n=numChildren)
        my.coverageoverSampbymNoUrb = coverage(thisTruth, designResNoUrb$overSampDat$Q10[,i],designResNoUrb$overSampDat$Q90[,i], doLogit=useLogit)
        my.lengthoverSampbymNoUrb = intervalWidth(designResNoUrb$overSampDat$Q10[,i], designResNoUrb$overSampDat$Q90[,i], logit=useLogit)
      }
      
      if("bym" %in% models) {
        my.biasoverSampbym = bias(thisTruth, designRes$overSampDat$mean[,i], logit=useLogit, my.var=(designRes$overSampDat$stddev[,i])^2)
        my.mseoverSampbym = mse(thisTruth, designRes$overSampDat$mean[,i], logit=useLogit, my.var=(designRes$overSampDat$stddev[,i])^2)
        my.dssoverSampbym = dss(thisTruth, designRes$overSampDat$mean[,i], (designRes$overSampDat$stddev[,i])^2)
        my.crpsoverSampbym = crpsNormal(thisTruth, designRes$overSampDat$mean[,i], (designRes$overSampDat$stddev[,i])^2, logit=useLogit, n=numChildren)
        my.coverageoverSampbym = coverage(thisTruth, designRes$overSampDat$Q10[,i],designRes$overSampDat$Q90[,i], doLogit=useLogit)
        my.lengthoverSampbym = intervalWidth(designRes$overSampDat$Q10[,i], designRes$overSampDat$Q90[,i], logit=useLogit)
      }
      
      if("spdeNoUrb" %in% models) {
        my.biasoverSampspdeNoUrb = bias(thisTruth, allresoverSamp$logit.est.spdeNoUrb, logit=useLogit, my.var=allresoverSamp$var.est.spdeNoUrb)
        my.mseoverSampspdeNoUrb = mse(thisTruth, allresoverSamp$logit.est.spdeNoUrb, logit=useLogit, my.var=allresoverSamp$var.est.spdeNoUrb)
        my.dssoverSampspdeNoUrb = dss(thisTruth, allresoverSamp$logit.est.spdeNoUrb, allresoverSamp$var.est.spdeNoUrb)
        my.crpsoverSampspdeNoUrb = crpsNormal(thisTruth, allresoverSamp$logit.est.spdeNoUrb, allresoverSamp$var.est.spdeNoUrb, logit=useLogit, n=numChildren)
        my.coverageoverSampspdeNoUrb = coverage(thisTruth, allresoverSamp$lower.spdeNoUrb, allresoverSamp$upper.spdeNoUrb, doLogit=useLogit)
        my.lengthoverSampspdeNoUrb = intervalWidth(allresoverSamp$lower.spdeNoUrb, allresoverSamp$upper.spdeNoUrb, logit=useLogit)
      }
      
      if("spde" %in% models) {
        my.biasoverSampspde = bias(thisTruth, allresoverSamp$logit.est.spde, logit=useLogit, my.var=allresoverSamp$var.est.spde)
        my.mseoverSampspde = mse(thisTruth, allresoverSamp$logit.est.spde, logit=useLogit, my.var=allresoverSamp$var.est.spde)
        my.dssoverSampspde = dss(thisTruth, allresoverSamp$logit.est.spde, allresoverSamp$var.est.spde)
        my.crpsoverSampspde = crpsNormal(thisTruth, allresoverSamp$logit.est.spde, allresoverSamp$var.est.spde, logit=useLogit, n=numChildren)
        my.coverageoverSampspde = coverage(thisTruth, allresoverSamp$lower.spde, allresoverSamp$upper.spde, doLogit=useLogit)
        my.lengthoverSampspde = intervalWidth(allresoverSamp$lower.spde, allresoverSamp$upper.spde, logit=useLogit)
      }
      
      if("direct" %in% models) {
        scoresDirectoverSamp <- rbind(scoresDirectoverSamp, 
                                      data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                                 bias=my.biasoverSampdirect, 
                                                 mse=my.mseoverSampdirect,
                                                 dss=my.dssoverSampdirect,
                                                 coverage=my.coverageoverSampdirect, 
                                                 var=mean(allresoverSamp$var.estdirect),
                                                 crps=my.crpsoverSampdirect, 
                                                 length=my.lengthoverSampdirect))
      }
      
      if("naive" %in% models) {
        scoresNaiveoverSamp<- rbind(scoresNaiveoverSamp, 
                                    data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                               bias=my.biasoverSampnaive, 
                                               mse=my.mseoverSampnaive,
                                               dss=my.dssoverSampnaive,
                                               coverage=my.coverageoverSampnaive, 
                                               var=mean(allresoverSamp$var.est),
                                               crps=my.crpsoverSampnaive, 
                                               length=my.lengthoverSampnaive))
      }
      
      if("mercer" %in% models) {
        scoresMerceroverSamp<- rbind(scoresMerceroverSamp, 
                                     data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                                bias=my.biasoverSampmercer, 
                                                mse=my.mseoverSampmercer,
                                                dss=my.dssoverSampmercer,
                                                coverage=my.coverageoverSampmercer, 
                                                var=mean(allresoverSamp$var.est.mercer),
                                                crps=my.crpsoverSampmercer, 
                                                length=my.lengthoverSampmercer))
      }
      
      if("bymNoUrb" %in% models) {
        scoresBYMNoUrboverSamp <- rbind(scoresBYMNoUrboverSamp, 
                                        data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                                   bias=my.biasoverSampbymNoUrb, 
                                                   mse=my.mseoverSampbymNoUrb,
                                                   dss=my.dssoverSampbymNoUrb,
                                                   coverage=my.coverageoverSampbymNoUrb, 
                                                   var=mean((designResNoUrb$overSampDat$stddev[,i])^2),
                                                   crps=my.crpsoverSampbymNoUrb, 
                                                   length=my.lengthoverSampbymNoUrb))
      }
      
      if("bym" %in% models) {
        scoresBYMoverSamp <- rbind(scoresBYMoverSamp, 
                                   data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                              bias=my.biasoverSampbym, 
                                              mse=my.mseoverSampbym,
                                              dss=my.dssoverSampbym,
                                              coverage=my.coverageoverSampbym, 
                                              var=mean((designRes$overSampDat$stddev[,i])^2),
                                              crps=my.crpsoverSampbym, 
                                              length=my.lengthoverSampbym))
      }
      
      if("spdeNoUrb" %in% models) {
        scoresSPDENoUrboverSamp <- rbind(scoresSPDENoUrboverSamp, 
                                         data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                                    bias=my.biasoverSampspdeNoUrb, 
                                                    mse=my.mseoverSampspdeNoUrb,
                                                    dss=my.dssoverSampspdeNoUrb,
                                                    coverage=my.coverageoverSampspdeNoUrb, 
                                                    var=mean(allresoverSamp$var.est.spdeNoUrb),
                                                    crps=my.crpsoverSampspdeNoUrb, 
                                                    length=my.lengthoverSampspdeNoUrb))
      }
      
      
      if("spde" %in% models) {
        scoresSPDEoverSamp <- rbind(scoresSPDEoverSamp, 
                                    data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                               bias=my.biasoverSampspde, 
                                               mse=my.mseoverSampspde,
                                               dss=my.dssoverSampspde,
                                               coverage=my.coverageoverSampspde, 
                                               var=mean(allresoverSamp$var.est.spde),
                                               crps=my.crpsoverSampspde, 
                                               length=my.lengthoverSampspde))
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
  
  # final table: 
  if("naive" %in% models)
    naive = apply(scoresNaiveSRS[, c("bias", "mse", "dss", "crps", "var", "coverage", "length")], 2, mean)
  if("direct" %in% models)
    direct = apply(scoresDirectSRS[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("mercer" %in% models)
    mercer = apply(scoresMercerSRS[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("bymNoUrb" %in% models)
    bymNoUrb = apply(scoresBYMNoUrbSRS[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("bym" %in% models)
    bym = apply(scoresBYMSRS[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("spdeNoUrb" %in% models)
    spdeNoUrb = apply(scoresSPDENoUrbSRS[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("spde" %in% models)
    spde = apply(scoresSPDESRS[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  idx = 1:9
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
  if("bymNoClust" %in% models)
    tab = rbind(tab, c(bymNoClust[idx]))
  if("bymNoUrbClust" %in% models)
    tab = rbind(tab, c(bymNoUrbClust[idx]))
  if("spde" %in% models)
    tab = rbind(tab, c(spde[idx]))
  if("spdeNoUrb" %in% models)
    tab = rbind(tab, c(spdeNoUrb[idx]))
  
  rownames(tab) = allNames[modelsI]
  
  print(xtable(tab, digits=3), 
        only.contents=TRUE, 
        include.colnames=TRUE,
        hline.after=NULL)
  
  if(produceFigures) {
    # compare all five
    pdf(paste0("figures/biasbyregion", sampling, ".pdf"), width=20, height=12)
    par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
    boxplot(bias~region, data=scoresDirect, at=seq(-1, 275, by=6), 
            col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
    boxplot(bias~region, data=scoresNaive, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
    boxplot(bias~region, data=scoresMercer, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
    boxplot(bias~region, data=scoresBYM, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
    boxplot(bias~region, data=scoresSPDE, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
    # axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirect$region[1:47])
    # axis(1, at=seq(0.5, 276.5, by=6), labels=scoresDirect$region[1:47])
    legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
           fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
    abline(h=0, lwd=2, col=2)
    dev.off()
    
    pdf(paste0("figures/crpsbyregion", sampling, ".pdf"), width=20, height=12)
    par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
    boxplot(crps~region, data=scoresDirect, at=seq(-1, 275, by=6), 
            col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
    boxplot(crps~region, data=scoresNaive, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
    boxplot(crps~region, data=scoresMercer, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
    boxplot(crps~region, data=scoresBYM, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
    boxplot(crps~region, data=scoresSPDE, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
    # axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirect$region[1:47])
    legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
           fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
    dev.off()
  }
  
  tab
}

getTruth = function(resultType = c("county", "region", "EA", "pixel"), eaDat) {
  resultType = match.arg(resultType)
  
  if(resultType == "county") {
    # compute truth based on superpopulation
    regions = sort(unique(eaDat$admin1))
    truthbycounty <- rep(NA, 47)
    childrenPerCounty = truthbycounty
    
    for(i in 1:47){
      super = eaDat[eaDat$admin1 == regions[i],]
      childrenPerCounty[i] = sum(super$numChildren)
      truthbycounty[i] <- sum(super$died)/childrenPerCounty[i]
    }
    truth = data.frame(admin1=regions, truth=truthbycounty, numChildren=childrenPerCounty)
  } else if(resultType == "pixel") {
    eaToPixel = eaDat$pixelI
    pixelsWithData = unique(eaToPixel)
    counties = sort(unique(eaDat$admin1))
    childrenPerPixel = aggregate(eaDat$numChildren, list(pixel=eaDat$pixelI), sum)
    urbanPixel = aggregate(eaDat$urban, list(pixel=eaDat$pixelI), function(x) {mean(x[1])})
    deathsPerPixel = aggregate(eaDat$died, list(pixel=eaDat$pixelI), sum)
    regions = childrenPerPixel$pixel # these are the pixels with enumeration areas in them
    # sort results by pixels with data
    sortI=match(pixelsWithData, regions)
    
    pixelToAdmin = match(popGrid$admin1[as.numeric(regions)], counties)
    
    truth = data.frame(pixel=regions, truth=deathsPerPixel$x / childrenPerPixel$x, countyI=pixelToAdmin, urban=urbanPixel$x, numChildren=childrenPerPixel$x)
    truth = truth[sortI,]
  } else if(resultType == "EA") {
    truth = data.frame(EA = 1:nrow(eaDat), truth = eaDat$died/eaDat$numChildren, urban=eaDat$urban, numChildren=eaDat$numChildren)
  } else if(resultType == "region") {
    regions = sort(unique(eaDat$region))
    truthbyregion <- rep(NA, 8)
    childrenPerRegion = truthbyregion
    
    for(i in 1:8){
      super = eaDat[eaDat$region == regions[i],]
      childrenPerRegion[i] = sum(super$numChildren)
      truthbyregion[i] <- sum(super$died)/childrenPerRegion[i]
    }
    truth = data.frame(admin1=regions, truth=truthbyregion, numChildren=childrenPerRegion)
  }
  
  truth
}