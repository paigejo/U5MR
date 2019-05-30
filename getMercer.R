# load direct estimates
# load("resultsDirectNaive.RData")

# source Mercer et al. code
source("mercer.R")
library(INLA)

# load a different 1 of these depending on whether a cluster effect should be included 
# in the simulation of the data or not (tausq is the cluster effect variance)
getMercer = function(tausq=.1^2, test=FALSE, margVar=0.15^2, gamma=-1) {
  if(!test)
    load(paste0("resultsDirectNaiveBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
           "HHoldVar0urbanOverSamplefrac0.RData"))
  else
    load(paste0("resultsDirectNaiveBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
           "HHoldVar0urbanOverSamplefrac0Test.RData"))
  
  # number of simulation scenarios
  n = 100
  mercerSRS = merceroverSamp = list()
  mercerSRSParEst = matrix(nrow=5, ncol=n)
  mercerSRSParSD = matrix(nrow=5, ncol=n)
  mercerSRSPar10 = matrix(nrow=5, ncol=n)
  mercerSRSPar50 = matrix(nrow=5, ncol=n)
  mercerSRSPar90 = matrix(nrow=5, ncol=n)
  merceroverSampParEst = matrix(nrow=5, ncol=n)
  merceroverSampParSD = matrix(nrow=5, ncol=n)
  merceroverSampPar10 = matrix(nrow=5, ncol=n)
  merceroverSampPar50 = matrix(nrow=5, ncol=n)
  merceroverSampPar90 = matrix(nrow=5, ncol=n)
  parNames = c("Intercept", "BYM2 Phi", "BYM2 Tot. Var", "BYM2 Spatial Var", "BYM2 iid Var")
  rownames(mercerSRSParEst) = parNames
  rownames(mercerSRSParSD) = parNames
  rownames(mercerSRSPar10) = parNames
  rownames(mercerSRSPar50) = parNames
  rownames(mercerSRSPar90) = parNames
  rownames(merceroverSampParEst) = parNames
  rownames(merceroverSampParSD) = parNames
  rownames(merceroverSampPar10) = parNames
  rownames(merceroverSampPar50) = parNames
  rownames(merceroverSampPar90) = parNames
  
  for(i in 1:n){
    print(i)
    
    # get predictive summary statistics
    tmpoverSamp = mercer_u1m2(directEstoverSamp[[i]]$logit.est, directEstoverSamp[[i]]$var.est, 
                             graph.path = "Kenyaadm1.graph")
    resoverSamp=   data.frame(admin1=directEstoverSamp[[i]]$admin1, 
                              u1m.mercer=expit(tmpoverSamp$summary.linear.predictor$mean),
                              lower.mercer=tmpoverSamp$summary.linear.predictor$"0.1quant",
                              upper.mercer=tmpoverSamp$summary.linear.predictor$"0.9quant",
                              logit.est.mercer=tmpoverSamp$summary.linear.predictor$mean, 
                              var.est.mercer=(tmpoverSamp$summary.linear.predictor$sd)^2)
    merceroverSamp[[i]] = resoverSamp
    
    
    tmpSRS = mercer_u1m2(directEstSRS[[i]]$logit.est, directEstSRS[[i]]$var.est, 
                        graph.path = "Kenyaadm1.graph")
    
    resSRS = data.frame(admin1=directEstSRS[[i]]$admin1,
                        u1m.mercer=expit(tmpSRS$summary.linear.predictor$mean),
                        lower.mercer=tmpSRS$summary.linear.predictor$"0.1quant",
                        upper.mercer=tmpSRS$summary.linear.predictor$"0.9quant",
                        logit.est.mercer=tmpSRS$summary.linear.predictor$mean, 
                        var.est.mercer=(tmpSRS$summary.linear.predictor$sd)^2)
    mercerSRS[[i]] = resSRS
    
    ## include parameter estimates in the tables (SRS)
    # intercept
    mercerSRSParEst[1, i] = tmpSRS$summary.fixed[,1]
    mercerSRSParSD[1, i] = tmpSRS$summary.fixed[,2]
    mercerSRSPar10[1, i] = tmpSRS$summary.fixed[,3]
    mercerSRSPar50[1, i] = tmpSRS$summary.fixed[,4]
    mercerSRSPar90[1, i] = tmpSRS$summary.fixed[,5]
    
    # BYM2 hyperparameter phi
    mercerSRSParEst[2, i] = tmpSRS$summary.hyperpar[2,1]
    mercerSRSParSD[2, i] = tmpSRS$summary.hyperpar[2,2]
    mercerSRSPar10[2, i] = tmpSRS$summary.hyperpar[2,3]
    mercerSRSPar50[2, i] = tmpSRS$summary.hyperpar[2,4]
    mercerSRSPar90[2, i] = tmpSRS$summary.hyperpar[2,5]
    
    ## transformed hyperparameters
    # sample the hyperparameters, using the marginals to improve the sampling
    out = inla.hyperpar.sample(1000, tmpSRS, improve.marginals=TRUE)
    transformFunction = function(x) {c(1/x[1], 1/x[1]*x[2], 1/x[1]*(1-x[2]))}
    transformedOut = apply(out, 1, transformFunction)
    
    # now calculate the summary statistics of the transformed BYM2 hyperparameters
    mercerSRSParEst[3:5, i] = rowMeans(transformedOut[1:3,])
    mercerSRSParSD[3:5, i] = apply(transformedOut[1:3,], 1, sd)
    mercerSRSPar10[3:5, i] = apply(transformedOut[1:3,], 1, quantile, probs=.1)
    mercerSRSPar50[3:5, i] = apply(transformedOut[1:3,], 1, quantile, probs=.5)
    mercerSRSPar90[3:5, i] = apply(transformedOut[1:3,], 1, quantile, probs=.9)
    
    ## include parameter estimates in the tables (urban oversampled)
    # intercept
    merceroverSampParEst[1, i] = tmpoverSamp$summary.fixed[,1]
    merceroverSampParSD[1, i] = tmpoverSamp$summary.fixed[,2]
    merceroverSampPar10[1, i] = tmpoverSamp$summary.fixed[,3]
    merceroverSampPar50[1, i] = tmpoverSamp$summary.fixed[,4]
    merceroverSampPar90[1, i] = tmpoverSamp$summary.fixed[,5]
    
    # BYM2 hyperparameter phi
    merceroverSampParEst[2, i] = tmpoverSamp$summary.hyperpar[2,1]
    merceroverSampParSD[2, i] = tmpoverSamp$summary.hyperpar[2,2]
    merceroverSampPar10[2, i] = tmpoverSamp$summary.hyperpar[2,3]
    merceroverSampPar50[2, i] = tmpoverSamp$summary.hyperpar[2,4]
    merceroverSampPar90[2, i] = tmpoverSamp$summary.hyperpar[2,5]
    
    ## transformed hyperparameters
    # sample the hyperparameters, using the marginals to improve the sampling
    out = inla.hyperpar.sample(1000, tmpoverSamp, improve.marginals=TRUE)
    transformFunction = function(x) {c(1/x[1], 1/x[1]*x[2], 1/x[1]*(1-x[2]))}
    transformedOut = apply(out, 1, transformFunction)
    
    # now calculate the summary statistics of the transformed BYM2 hyperparameters
    merceroverSampParEst[3:5, i] = rowMeans(transformedOut[1:3,])
    merceroverSampParSD[3:5, i] = apply(transformedOut[1:3,], 1, sd)
    merceroverSampPar10[3:5, i] = apply(transformedOut[1:3,], 1, quantile, probs=.1)
    merceroverSampPar50[3:5, i] = apply(transformedOut[1:3,], 1, quantile, probs=.5)
    merceroverSampPar90[3:5, i] = apply(transformedOut[1:3,], 1, quantile, probs=.9)
  }
  
  # generate summary statistics about the parameters
  mercerSRSPar = data.frame(list(Est=rowMeans(mercerSRSParEst), 
                                 SD=rowMeans(mercerSRSParSD), 
                                 Q10=rowMeans(mercerSRSPar10), 
                                 Q50=rowMeans(mercerSRSPar50), 
                                 Q90=rowMeans(mercerSRSPar90)))
  merceroverSampPar = data.frame(list(Est=rowMeans(merceroverSampParEst), 
                                 SD=rowMeans(merceroverSampParSD), 
                                 Q10=rowMeans(merceroverSampPar10), 
                                 Q50=rowMeans(merceroverSampPar50), 
                                 Q90=rowMeans(merceroverSampPar90)))
  
  if(!test)
    save(merceroverSamp, mercerSRS, mercerSRSPar, merceroverSampPar, file=paste0("resultsMercerBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                                                "HHoldVar0urbanOverSamplefrac0.RData"))
  else
    save(merceroverSamp, mercerSRS, mercerSRSPar, merceroverSampPar, file=paste0("resultsMercerBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                                                "HHoldVar0urbanOverSamplefrac0Test.RData"))
  
  invisible(NULL)
}

# leave one county of data out at a time in order to validate the mercer model versus county level direct estimates
validateMercer = function(dat=ed, directLogitEsts, directLogitVars, directVars, 
                          counties=sort(unique(poppc$County))) {
  
  # fit the full model once, calculating certain validation scores
  print("Fitting full model")
  modelFit = mercer_u1m2(directLogitEsts, directLogitVars, 
                    graph.path = "Kenyaadm1.graph", doValidation = TRUE)
  cpo = modelFit$cpo$cpo
  dic = modelFit$dic$dic
  waic = modelFit$waic$waic
  
  # now do leave one county out across validation
  for(i in 1:length(counties)) {
    if(i %% 10 == 1)
      print(paste0("Fitting model with data from county ", i, "/", length(counties), " left out"))
    thisCountyName = counties[i]
    
    # fit model, get all predictions for each areal level and each posterior sample
    fit = mercer_u1m2(directLogitEsts, directLogitVars, 
                      graph.path = "Kenyaadm1.graph", 
                      previousResult=modelFit, predCountyI=i)
    
    # get predictive distribution for the left out county
    res = data.frame(admin1=counties,
                     est.mercer=expit(fit$summary.linear.predictor$mean),
                     lower.mercer=fit$summary.linear.predictor$"0.1quant",
                     upper.mercer=fit$summary.linear.predictor$"0.9quant",
                     logit.est.mercer=fit$summary.linear.predictor$mean, 
                     var.est.mercer=(fit$summary.linear.predictor$sd)^2)
    
    thisCountyPreds = res[i,]
    
    if(i == 1) {
      countyPreds = thisCountyPreds
    } else {
      countyPreds = rbind(countyPreds, thisCountyPreds)
    }
  }
  
  # calculate weights further predictions based on inverse direct estimate variances, calculate validation scores
  weights = 1 / directVars
  weights = weights / sum(weights)
  logitWeights = 1 / directLogitVars
  logitWeights = logitWeights / sum(logitWeights)
  theseScores = getValidationScores(directLogitEsts, directLogitVars, 
                                    countyPreds$logit.est.mercer, countyPreds$var.est.mercer, 
                                    weights=weights, logitWeights=logitWeights, usePearson=TRUE)
  
  # now calculate scoring rules on the cluster level
  countyI = match(dat$admin1, counties)
  theseScoresCluster = getValidationScores(dat$y / dat$n, rep(0, nrow(dat)), directEsts = dat$y / dat$n, 
                                           countyPreds$logit.est.mercer[countyI], countyPreds$var.est.mercer[countyI], 
                                           usePearson=TRUE)
  
  # compile scores c("MSE", "CPO", "CRPS", "logScore")
  scores = theseScores$scores
  pit = theseScores$pit
  scores = data.frame(MSE=scores$MSE, Biassq=scores$`Bias^2`, Var=scores$Var, CPO=scores$CPO, CRPS=scores$CRPS, logScore=scores$logScore)
  
  scoresCluster = theseScoresCluster$scores
  pitCluster = theseScoresCluster$pit
  scoresCluster = data.frame(DIC=dic, WAIC=waic, MSE=scoresCluster$MSE, Biassq=scoresCluster$`Bias^2`, Var=scoresCluster$Var, CPO=mean(cpo, na.rm=TRUE), CRPS=scoresCluster$CRPS, logScore=scoresCluster$logScore)
  
  list(scores=scores, pit=pit, pitCluster=pitCluster, scoresCluster=scoresCluster)
}

##### the below code appeared to be a copy of the above code, so I commented it out:
# # load a different 1 of these depending on whether a cluster effect should be included 
# # in the simulation of the data or not (tausq is the cluster effect variance)
# getMercer = function(tausq=.1^2, test=FALSE) {
#   # load("resultsDirectNaiveTausq0.RData")
#   if(!test)
#     load(paste0("resultsDirectNaiveTausq", round(tausq, 4), ".RData"))
#   else
#     load(paste0("resultsDirectNaiveTausq", round(tausq, 4), "Test.RData"))
#   
#   # number of simulation scenarios
#   n = 100
#   mercerSRS = merceroverSamp = list()
#   
#   for(i in 1:100){
#     print(i)
#     tmpoverSamp = mercer_u1m(directEstoverSamp[[i]]$logit.est, directEstoverSamp[[i]]$var.est, 
#                              graph.path = "Kenyaadm1.graph")
#     resoverSamp=   data.frame(admin1=directEstoverSamp[[i]]$admin1, 
#                               u1m.mercer=expit(tmpoverSamp$summary.linear.predictor$mean),
#                               lower.mercer=tmpoverSamp$summary.linear.predictor$"0.1quant",
#                               upper.mercer=tmpoverSamp$summary.linear.predictor$"0.9quant",
#                               logit.est.mercer=tmpoverSamp$summary.linear.predictor$mean, 
#                               var.est.mercer=(tmpoverSamp$summary.linear.predictor$sd)^2)
#     merceroverSamp[[i]] = resoverSamp
#     
#     
#     tmpSRS = mercer_u1m(directEstSRS[[i]]$logit.est, directEstSRS[[i]]$var.est, 
#                         graph.path = "Kenyaadm1.graph")
#     
#     resSRS = data.frame(admin1=directEstSRS[[i]]$admin1,
#                         u1m.mercer=expit(tmpSRS$summary.linear.predictor$mean),
#                         lower.mercer=tmpSRS$summary.linear.predictor$"0.1quant",
#                         upper.mercer=tmpSRS$summary.linear.predictor$"0.9quant",
#                         logit.est.mercer=tmpSRS$summary.linear.predictor$mean, 
#                         var.est.mercer=(tmpSRS$summary.linear.predictor$sd)^2)
#     mercerSRS[[i]] = resSRS
#   }
#   
#   # save(merceroverSamp, mercerSRS, file="resultsMercer.RData")
#   # save(merceroverSamp, mercerSRS, file="resultsMercerTausq0.RData")
#   if(!test)
#     save(merceroverSamp, mercerSRS, file=paste0("resultsMercerTausq", round(tausq, 4), ".RData"))
#   else
#     save(merceroverSamp, mercerSRS, file=paste0("resultsMercerTausq", round(tausq, 4), "test.RData"))
#   
#   invisible(NULL)
# }