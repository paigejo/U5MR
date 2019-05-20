# validate the smoothing models by leaving out data from one county at a time

validateExampleHelper = function(dat=ed, resultNameRoot="Ed", leftOutCountyI, counties=sort(unique(poppc$County))) {
  resultNameRootLower = tolower(resultNameRoot)
  
  # remove data from the designated county
  leftOutCounty = counties[leftOutCountyI]
  dat = dat[dat$admin1 != leftOutCounty,]
  
  # no need to run the direct and naive examples, since they cannot make predictions for a county that does not have data
  
  ##### run Mercer et al. model
  source("mercer.R")
  tmpResults = mercer_u1m2(directEstResults$logit.est, directEstResults$var.est, 
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
  runBYM2Dat(dat, includeUrbanRural = FALSE, includeCluster = FALSE, fileNameRoot=resultNameRoot)
  runBYM2Dat(dat, includeUrbanRural = FALSE, includeCluster = TRUE, fileNameRoot=resultNameRoot)
  runBYM2Dat(dat, includeUrbanRural = TRUE, includeCluster = FALSE, fileNameRoot=resultNameRoot)
  runBYM2Dat(dat, includeUrbanRural = TRUE, includeCluster = TRUE, fileNameRoot=resultNameRoot)
  
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