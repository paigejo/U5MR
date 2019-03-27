# function for analyzing example datasets

# first name elements of ed to be the same as the corresponding elements of the simulated datasets

generateExampleResults = function(dat=ed, resultNameRoot="Ed") {
  resultNameRootLower = tolower(resultNameRoot)
  
  ##### first generate results for direct and naive models
  startIDat = 1
  for(i in 1:nrow(dat)) {
    if(i %% 100 == 1)
      print(paste0("i: ", i))
    
    tmpDat = extendDataDat(dat[i,], v001 = i)
    
    # initialize the data frames on the first iteration only
    if(i == 1) {
      resDat = as.data.frame(matrix(nrow=sum(dat$n), ncol=ncol(tmpDat) + 1))
    }
    
    # append to the data frames
    names(resDat) = c(names(tmpDat), "regionRural")
    endIDat = startIDat + nrow(tmpDat) - 1
    resDat[startIDat:endIDat, 1:ncol(tmpDat)] = tmpDat
    
    # update row index
    startIDat = endIDat + 1
  }
  
  # add in RegionRural interaction
  resDat$regionRural <- with(resDat, interaction(admin1, urban), drop=TRUE)
  
  # save the resulting data frame
  save(resDat, file=paste0("data4direct", resultNameRoot, ".RData"))
  
  # start analysis using naive approach and computing direct estimates
  directEstDat = naiveDat = list()
  
  # analyse the unstratified sampling scenario
  dat_obj2 = resDat
  res2 = defineSurveyDat(dat_obj2, 
                         stratVar=dat_obj2$regionRural,
                         useSamplingWeights = TRUE)
  directEstResults = res2
  defineSurvey
  resnA2 = run_naiveDat(dat_obj2)
  
  naiveResults = resnA2
  
  # save results
  save(directEstResults, naiveResults, file=paste0("resultsDirectNaive", resultNameRoot, ".RData"))
  
  ##### run Mercer et al. model
  source("mercer.R")
  tmpResults = mercer_u1m(directEstResults$logit.est, directEstResults$var.est, 
                          graph.path = "Kenyaadm1.graph")
  
  res = data.frame(admin1=directEstResults$admin1,
                   est.mercer=expit(tmpResults$summary.linear.predictor$mean),
                   lower.mercer=tmpResults$summary.linear.predictor$"0.1quant",
                   upper.mercer=tmpResults$summary.linear.predictor$"0.9quant",
                   logit.est.mercer=tmpResults$summary.linear.predictor$mean, 
                   var.est.mercer=(tmpResults$summary.linear.predictor$sd)^2)
  mercerResults = res
  
  save(mercerResults, file=paste0("resultsMercer", resultNameRoot, ".RData"))
  
  ##### run BYM models
  source("designBased.R")
  runBYM2Dat(dat, includeUrbanRural = FALSE, includeCluster = FALSE)
  runBYM2Dat(dat, includeUrbanRural = FALSE, includeCluster = TRUE)
  runBYM2Dat(dat, includeUrbanRural = TRUE, includeCluster = FALSE)
  runBYM2Dat(dat, includeUrbanRural = TRUE, includeCluster = TRUE)
  
  ##### run SPDE 
  argList = list(list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  for(i in 1:length(argList)) {
    args = argList[[i]]
    spdeResults = do.call("resultsSPDEDat", args)
    includeClustEffect = args$includeClustEffect
    urbanEffect = args$urbanEffect
    fileName = paste0("resultsSPDE", resultNameRootLower, "_includeClustEffect", includeClustEffect, 
                      "_urbanEffect", urbanEffect, ".RData")
    save(spdeResults, file=fileName)
  }
  invisible(NULL)
}