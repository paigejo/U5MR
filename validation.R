# validate the smoothing models by leaving out data from one county at a time

validateExample = function(dat=ed, resultNameRoot="Ed", directEstResults=NULL, counties=sort(unique(poppc$County)), 
                           startFrom=0) {
  
  if(is.null(directEstResults)) {
    print("Loading Direct estimates results...")
    out = load(paste0("resultsDirectNaive", resultNameRoot, ".RData"))
    # directEstResults should now be loaded into the environment
  }
  
  # no need to run the direct and naive examples, since they cannot make predictions for a county that does not have data
  
  # calculate the variance of the direct estimates assuming Gaussianity on the logit scale to get the validation weights
  directMoments = logitNormMoments(cbind(directEstResults$logit.est, directEstResults$var.est))
  directEstResults$varProbScale = directMoments[,2]
  
  ##### run Mercer et al. model
  if(startFrom <= 0) {
    print("Validating Mercer model...")
    source("mercer.R")
    mercerResults = validateMercer(dat=dat, directEstResults$logit.est, directEstResults$var.est, directEstResults$varProbScale)
    
    save(mercerResults, file=paste0("resultsMercer", resultNameRoot, "ValidationAll.RData"))
  }
  else {
    print("Loading Smoothed Direct results...")
    load(file=paste0("resultsMercer", resultNameRoot, "ValidationAll.RData"))
  }
  
  ##### run BYM models
  source("designBased.R")
  if(startFrom <= 1) {
    print("Validating BYM2 I model...")
    bym2I = validateBYM2Dat(directEstResults$logit.est, directEstResults$var.est, directEstResults$varProbScale, 
                            dat, saveResults=TRUE, includeUrbanRural = FALSE, includeCluster = FALSE)$bym2Results
  }
  else {
    print("Loading BYM2 I results...")
    load(paste0('bym2', resultNameRoot, 'ValidationAllUrbRur', FALSE, 'Cluster', FALSE, '.RData'))
    bym2I = bym2Results
  }
  
  if(startFrom <= 2) {
    print("Validating BYM2 II model...")
    out = validateBYM2Dat(directEstResults$logit.est, directEstResults$var.est, directEstResults$varProbScale, 
                             dat, saveResults=TRUE, includeUrbanRural = FALSE, includeCluster = TRUE)
    bym2II = out$bym2Results
    bym2IIMod = out$bym2ResultsMod
  }
  else {
    print("Loading BYM2 II results...")
    load(paste0('bym2', resultNameRoot, 'ValidationAllUrbRur', FALSE, 'Cluster', TRUE, '.RData'))
    bym2II = bym2Results
    bym2IIMod = bym2ResultsMod
  }
  
  if(startFrom <= 3) {
    print("Validating BYM2 III model...")
    bym2III = validateBYM2Dat(directEstResults$logit.est, directEstResults$var.est, directEstResults$varProbScale, 
                              dat, saveResults=TRUE, includeUrbanRural = TRUE, includeCluster = FALSE)
  }
  else {
    print("Loading BYM2 III results...")
    load(paste0('bym2', resultNameRoot, 'ValidationAllUrbRur', TRUE, 'Cluster', FALSE, '.RData'))
    bym2III = bym2Results
  }
  
  if(startFrom <= 4) {
    print("Validating BYM2 IV model...")
    out = validateBYM2Dat(directEstResults$logit.est, directEstResults$var.est, directEstResults$varProbScale, 
                             dat, saveResults=TRUE, includeUrbanRural = TRUE, includeCluster = TRUE)
    bym2IV = out$bym2Results
    bym2IVMod = out$bym2ResultsMod
  }
  else {
    print("Loading BYM2 IV results...")
    load(paste0('bym2', resultNameRoot, 'ValidationAllUrbRur', TRUE, 'Cluster', TRUE, '.RData'))
    bym2IV = bym2Results
    bym2IVMod = bym2ResultsMod
  }
  
  ##### run SPDE models
  if(startFrom <= 5) {
    print("Validating SPDE I model...")
    spdeI = validateSPDEDat(clustDat = dat, directLogitEsts=directEstResults$logit.est, directLogitVars=directEstResults$var.est, 
                            directVars=directEstResults$varProbScale, includeClustEffect = FALSE, urbanEffect = FALSE, saveResults = TRUE)
  }
  else {
    print("Loading SPDE I results...")
    fileName = paste0("resultsSPDE", resultNameRoot, "ValidationAll", "_includeClustEffect", FALSE, 
                      "_urbanEffect", FALSE, ".RData")
    load(fileName)
    spdeI = spdeResults
  }
  
  if(startFrom <= 6) {
    print("Validating SPDE II model...")
    spdeII = validateSPDEDat(clustDat = dat, directLogitEsts=directEstResults$logit.est, directLogitVars=directEstResults$var.est, 
                             directVars=directEstResults$varProbScale, includeClustEffect = FALSE, urbanEffect = TRUE, saveResults = TRUE)
  }
  else {
    print("Loading SPDE II results...")
    fileName = paste0("resultsSPDE", resultNameRoot, "ValidationAll", "_includeClustEffect", FALSE, 
                      "_urbanEffect", TRUE, ".RData")
    load(fileName)
    spdeII = spdeResults
  }
  
  if(startFrom <= 7) {
    print("Validating SPDE III model...")
    spdeIII = validateSPDEDat(clustDat = dat, directLogitEsts=directEstResults$logit.est, directLogitVars=directEstResults$var.est, 
                              directVars=directEstResults$varProbScale, includeClustEffect = TRUE, urbanEffect = FALSE, saveResults = TRUE)
  }
  else {
    print("Loading SPDE III results...")
    fileName = paste0("resultsSPDE", resultNameRoot, "ValidationAll", "_includeClustEffect", TRUE, 
                      "_urbanEffect", FALSE, ".RData")
    load(fileName)
    spdeIII = spdeResults
  }
  
  if(startFrom <= 8) {
    print("Validating SPDE IV model...")
    spdeIV = validateSPDEDat(clustDat = dat, directLogitEsts=directEstResults$logit.est, directLogitVars=directEstResults$var.est, 
                             directVars=directEstResults$varProbScale, includeClustEffect = TRUE, urbanEffect = TRUE, saveResults = TRUE)
  }
  else {
    print("Loading SPDE IV results...")
    fileName = paste0("resultsSPDE", resultNameRoot, "ValidationAll", "_includeClustEffect", TRUE, 
                      "_urbanEffect", TRUE, ".RData")
    load(fileName)
    spdeIV = spdeResults
  }
  
  ##### aggregate results
  print("Aggregating and saving results...")
  tab = rbind(mercerResults$scores, 
              bym2I$scores, 
              bym2II$scores, 
              bym2IIMod$scores, 
              bym2III$scores, 
              bym2IV$scores, 
              bym2IVMod$scores, 
              spdeI$scores, 
              spdeII$scores, 
              spdeIII$scores, 
              spdeIV$scores)
  rownames(tab) = c("Smoothed Direct", "BYM2 I", "BYM2 II", "BYM2 II'", "BYM2 III", "BYM2 IV", "BYM2 IV'", "SPDE I", "SPDE II", "SPDE III", "SPDE IV")
  PITs = list(mercer=mercerResults$pit, 
              bym2I=bym2I$pit, 
              bym2II=bym2II$pit, 
              bym2IIMod=bym2IIMod$pit, 
              bym2III=bym2III$pit, 
              bym2IV=bym2IV$pit, 
              bym2IVMod=bym2IVMod$pit, 
              spdeI=spdeI$pit, 
              spdeII=spdeII$pit, 
              spdeIII=spdeIII$pit, 
              spdeIV=spdeIV$pit)
  tabCluster = rbind(mercerResults$scoresCluster, 
              bym2I$scoresCluster, 
              bym2II$scoresCluster, 
              bym2IIMod$scoresCluster, 
              bym2III$scoresCluster, 
              bym2IV$scoresCluster, 
              bym2IVMod$scoresCluster, 
              spdeI$scoresCluster, 
              spdeII$scoresCluster, 
              spdeIII$scoresCluster, 
              spdeIV$scoresCluster)
  rownames(tabCluster) = c("Smoothed Direct", "BYM2 I", "BYM2 II", "BYM2 II'", "BYM2 III", "BYM2 IV", "BYM2 IV'", "SPDE I", "SPDE II", "SPDE III", "SPDE IV")
  PITsCluster = list(mercer=mercerResults$pitCluster, 
              bym2I=bym2I$pitCluster, 
              bym2II=bym2II$pitCluster, 
              bym2IIMod=bym2IIMod$pitCluster, 
              bym2III=bym2III$pitCluster, 
              bym2IV=bym2IV$pitCluster, 
              bym2IVMod=bym2IVMod$pitCluster, 
              spdeI=spdeI$pitCluster, 
              spdeII=spdeII$pitCluster, 
              spdeIII=spdeIII$pitCluster, 
              spdeIV=spdeIV$pitCluster)
  
  # save and return results
  validationResults = list(scores=tab, pit=PITs, scoresCluster=tabCluster, pitCluster=PITsCluster)
  save(validationResults, file=paste0("resultsValidationAll", resultNameRoot, ".RData"))
  invisible(validationResults)
}

# print out validation results and plot the PITs
printValidationResults = function(resultNameRoot="Ed") {
  # first load the validation results
  out = load(paste0("resultsValidationAll", resultNameRoot, ".RData"))
  modelNames = rownames(validationResults$scores)
  
  ## plot the PITs
  for(i in 1:length(modelNames)) {
    hist(validationResults$pit[[i]], main=paste0("Histogram of ", modelNames[i], " PITs"), 
         xlab="PIT", freq=FALSE, breaks=30, col="skyblue")
  }
  
  ## modify the BYM2 model names 
  tab = validationResults$scores
  
  ## print out the scoring rule table
  xtable(validationResults$scores, digits=c(1, 1, 1, rep(3, ncol(validationResults$scores)-2)), 
         display=c("s", "f", "f", rep("fg", ncol(validationResults$scores)-2)))
}





