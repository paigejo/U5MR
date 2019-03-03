# load("simDat.RData")
# library(profvis)
# library(logitnorm)
source("setup.R")
# setwd("~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/")

##### the code below does not use the same enumeration areas for each simulation, 
##### which was why was commented out


# 
# # number of datasets to be generated
# numData <- 100
# 
# set.seed(580252)
# my.seeds <- round(runif(numData)*100000)
# 
# dataSets = list()
# for(i in 1:numData){
#   print(i)
#   beta0 = -2
#   margVar = .15^2
#   tausq = .1^2
#   gamma = -.5
#   HHoldVar = 0
#   # HHoldVar = .3^2
#   tmp <- simDat(kenyaDat, beta0=beta0, margVar=margVar, tausq=tausq, gamma=gamma, seed=my.seeds[i])
#   # add the vector of sampling weights to the dataset 
#   #samplingWeight = 1/(table(tmp$clustDat$admin1)/table(tmp$eaDat$admin1))
#   #tmp$clustDat$samplingWeight = samplingWeight
#   dataSets[[i]] = tmp
# }
# 
# save(dataSets, file=paste0("simData4analysisBeta", round(beta0, 2), "margVar", round(margVar, 2), "tausq", 
#                            round(tausq, 2), "gamma", round(gamma, 2), "HHoldVar", HHoldVar, ".RData"))
# 

## CAUTION!!!!
#Warning messages:
# 1: In doTryCatch(return(expr), name, parentenv, handler) :
#  restarting interrupted promise evaluation
#2: In doTryCatch(return(expr), name, parentenv, handler) :
#  restarting interrupted promise evaluation



# unlike the above script, this script holds the EA data for each simulation fixed, 
# each simulation instead varying which clusters were sampled
# urbanOverSample: within any county, any individual urban EA is urbanOverSample times 
#                  as likely to be sampled than any individual rural EA to be a cluster
# nsim=100

# generate the empirical distributions and save them
# wd = getwd()
# setwd("~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/")
# empiricalDistributions = getSurveyEmpiricalDistributions2()
# save(empiricalDistributions, file="empiricalDistributions.RData")
# save(empiricalDistributions, file="~/git/U5MR/empiricalDistributions.RData")
# setwd(wd)


# simulate and save datasets used for the simulation study with the given model parameters
# NOTE: paired with the dataset using the passed parameters will be another dataset from the 
#       same model without a nugget/cluster effect
# nsim: number of surveys taken from the true latent population in the standard size survey collections
# nsimBig: number of surveys taken from the true latent population in the large size survey collections
# seeds: random number seeds used for making the latent population and generating surveys respectively
# beta0: latent gaussian model intercept
# margVar: marginal variance of the spatial field
# tausq: the nugget/cluster effect variance
# gamma: latent gaussian model urban effect
# HHoldVar: household effect variance
# effRange: spatial range
# urbanOverSamplefrac: the proportion with which to inflate the amount of urban samples in the surveys
generateSimDataSets = function(nsim=100, nsimBig = 250, seeds=c(580252, 1234), beta0 = -1.75, margVar = .15^2, 
                               tausq = .1^2, gamma = -1, HHoldVar = 0, effRange = 150, 
                               urbanOverSamplefrac = 0) {
  set.seed(seeds[1])
  wd = getwd()
  setwd("~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/")
  
  # make strings representing the simulation with different numbers of effects
  # full model
  dataID = paste0("Beta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                  round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, 
                  "urbanOverSamplefrac", round(urbanOverSamplefrac, 4))
  # no cluster effect
  dataID0 = paste0("Beta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                   round(0, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, 
                   "urbanOverSamplefrac", round(urbanOverSamplefrac, 4))
  # no cluster or urban effects (constant plus the spatial effect only)
  dataID02 = paste0("Beta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                   round(0, 4), "gamma", round(0, 4), "HHoldVar", HHoldVar, 
                   "urbanOverSamplefrac", round(urbanOverSamplefrac, 4))
  # constant effect only
  dataID03 = paste0("Beta", round(beta0, 4), "margVar", round(0, 4), "tausq", 
                   round(0, 4), "gamma", round(0, 4), "HHoldVar", HHoldVar, 
                   "urbanOverSamplefrac", round(urbanOverSamplefrac, 4))
  
  # there should be 1 true data set, but many simulated cluster samples
  load("empiricalDistributions.RData")
  simulatedEAs = simDatEmpirical(empiricalDistributions, kenyaEAs, clustDat=NULL, nsim=1, 
                                 beta0=beta0, margVar=margVar, urbanOverSamplefrac=urbanOverSamplefrac, 
                                 tausq=tausq, gamma=gamma, HHoldVar=HHoldVar, effRange=effRange)
  kenyaEAs = simulatedEAs$eaDat
  kenyaEAs$eaIs = 1:nrow(kenyaEAs)
  kenyaEAsLong = kenyaEAs[rep(1:nrow(kenyaEAs), kenyaEAs$nHH),]
  
  set.seed(seeds[2])
  # simulate the cluster sampling and add to the data sets
  overSampClustDat = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsimBig, NULL, urbanOverSamplefrac, verbose=FALSE)
  overSampClustDat = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsimBig, NULL, urbanOverSamplefrac, verbose=FALSE)
  clustList = genAndreaFormatFromEAIs(simulatedEAs$eaDat, overSampClustDat$eaIs, overSampClustDat$sampleWeights)
  overSampDat = list(eaDat=kenyaEAs, clustDat=clustList)
  
  overSampClustDatTest = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsimBig, NULL, urbanOverSamplefrac, fixedPerStrata=TRUE, nPerStrata=3, verbose=FALSE)
  clustListTest = genAndreaFormatFromEAIs(kenyaEAs, overSampClustDatTest$eaIs, overSampClustDatTest$sampleWeights)
  overSampDatTest = list(eaDat=kenyaEAs, clustDat=clustListTest)
  
  SRSClustDat = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsimBig, NULL, SRS=TRUE, verbose=FALSE)
  clustList = genAndreaFormatFromEAIs(kenyaEAs, SRSClustDat$eaIs, SRSClustDat$sampleWeights)
  SRSDat = list(eaDat=kenyaEAs, clustDat=clustList) # the only thing different is the sampling of the clusters
  
  SRSClustDatTest = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsimBig, NULL, fixedPerStrata=TRUE, nPerStrata=3, SRS=TRUE, verbose=FALSE)
  clustListTest = genAndreaFormatFromEAIs(kenyaEAs, SRSClustDatTest$eaIs, SRSClustDatTest$sampleWeights)
  SRSDatTest = list(eaDat=kenyaEAs, clustDat=clustListTest) # the only thing different is the sampling of the clusters
  
  # plot the first simulation of the over sampled and simple random sample data sets
  clustDat = SRSDat$clustDat[[1]]
  # clustDat = overSampDat$clustDat[[1]]
  eaDat = overSampDat$eaDat
  pdf(paste0("figures/exampleSRSSimulation", dataID, ".pdf"), width=8, height=8)
  par(mfrow =c(2, 2))
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = clustDat$numChildren
  obsCounts = clustDat$died
  zlim = c(0, quantile(c(eaDat$died/eaDat$numChildren, clustDat$died/clustDat$numChildren, 
                         eaDat$trueProbDeath), probs=.975))
  quilt.plot(eaDat$east, eaDat$north, eaDat$died/eaDat$numChildren, main="All Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, obsCounts/obsNs, main="Sample Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(eaDat$east, eaDat$north, eaDat$trueProbDeath, main="All True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, clustDat$trueProbDeath, main="Sample True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  dev.off()
  
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID, "Big.RData"))
  overSampDat = overSampDatTest
  SRSDat = SRSDatTest
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID, "TestBig.RData"))
  out = load(paste0("simDataMulti", dataID, "Big.RData"))
  
  # Now take only the first nsim simulations from the "big" dataset
  overSampDat$clustDat = overSampDat$clustDat[1:nsim]
  SRSDat$clustDat = SRSDat$clustDat[1:nsim]
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID, ".RData"))
  overSampDat = overSampDatTest
  SRSDat = SRSDatTest
  overSampDat$clustDat = overSampDat$clustDat[1:nsim]
  SRSDat$clustDat = SRSDat$clustDat[1:nsim]
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID, "Test.RData"))
  
  # reload the data
  out = load(paste0("simDataMulti", dataID, "Big.RData"))
  
  # Now simulate the data without a cluster effect but with the same underlying probability surface otherwise
  tausq = 0
  overSampDat$eaDat$trueProbDeath = overSampDat$eaDat$trueProbDeathNoNug
  SRSDat$eaDat$trueProbDeath = SRSDat$eaDat$trueProbDeathNoNug
  overSampDat$eaDat$died = rbinom(nrow(overSampDat$eaDat), overSampDat$eaDat$numChildren, overSampDat$eaDat$trueProbDeathNoNug)
  SRSDat$eaDat$died = overSampDat$eaDat$died
  
  overSampDatTest$eaDat$trueProbDeath = overSampDatTest$eaDat$trueProbDeathNoNug
  SRSDatTest$eaDat$trueProbDeath = SRSDatTest$eaDat$trueProbDeathNoNug
  overSampDatTest$eaDat$died = rbinom(nrow(overSampDatTest$eaDat), overSampDatTest$eaDat$numChildren, overSampDatTest$eaDat$trueProbDeathNoNug)
  SRSDatTest$eaDat$died = overSampDatTest$eaDat$died
  for(i in 1:nsimBig) {
    overSampDat$clustDat[[i]]$trueProbDeath = overSampDat$clustDat[[i]]$trueProbDeathNoNug
    SRSDat$clustDat[[i]]$trueProbDeath = SRSDat$clustDat[[i]]$trueProbDeathNoNug
    overSampDatTest$clustDat[[i]]$trueProbDeath = overSampDatTest$clustDat[[i]]$trueProbDeathNoNug
    SRSDatTest$clustDat[[i]]$trueProbDeath = SRSDatTest$clustDat[[i]]$trueProbDeathNoNug
    
    overSampDat$clustDat[[i]]$died = overSampDat$eaDat$died[overSampClustDat$eaIs[,i]]
    SRSDat$clustDat[[i]]$died = SRSDat$eaDat$died[SRSClustDat$eaIs[,i]]
    overSampDatTest$clustDat[[i]]$died = overSampDatTest$eaDat$died[overSampClustDatTest$eaIs[,i]]
    SRSDatTest$clustDat[[i]]$died = SRSDatTest$eaDat$died[SRSClustDatTest$eaIs[,i]]
  }
  
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID0, "Big.RData"))
  overSampDat = overSampDatTest
  SRSDat = SRSDatTest
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID0, "TestBig.RData"))
  load(paste0("simDataMulti", dataID0, "Big.RData"))
  
  # Again, take only the first nsim simulations from the "big"" dataset
  overSampDat$clustDat = overSampDat$clustDat[1:nsim]
  SRSDat$clustDat = SRSDat$clustDat[1:nsim]
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID0, ".RData"))
  overSampDat = overSampDatTest
  SRSDat = SRSDatTest
  overSampDat$clustDat = overSampDat$clustDat[1:nsim]
  SRSDat$clustDat = SRSDat$clustDat[1:nsim]
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID0, "Test.RData"))
  
  # clustDat = SRSDat$clustDat[[1]]
  # clustDat = SRSDatTest$clustDat[[1]]
  # clustDat = overSampDat$clustDat[[1]]
  clustDat = overSampDatTest$clustDat[[1]]
  eaDat = overSampDat$eaDat
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = clustDat$numChildren
  obsCounts = clustDat$died
  
  pdf(paste0("figures/exampleOverSampTestSimulationNoNug", dataID0, ".pdf"), width=8, height=8)
  par(mfrow =c(2, 2))
  zlim = c(0, quantile(c(eaDat$died/eaDat$numChildren, clustDat$died/clustDat$numChildren, 
                         eaDat$trueProbDeath), probs=.975))
  quilt.plot(eaDat$east, eaDat$north, eaDat$died/eaDat$numChildren, main="All Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, obsCounts/obsNs, main="Sample Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(eaDat$east, eaDat$north, eaDat$trueProbDeath, main="All True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, clustDat$trueProbDeath, main="Sample True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  dev.off()
  
  ##### do the same thing but removing the urban effect
  # reload the data
  out = load(paste0("simDataMulti", dataID, "Big.RData"))
  
  # Now simulate the data without a cluster effect or urban effect but with the same underlying probability surface otherwise
  tausq = 0
  overSampDat$eaDat$trueProbDeathNoNug = expit(logit(overSampDat$eaDat$trueProbDeathNoNug) - gamma*overSampDat$eaDat$urban)
  overSampDat$eaDat$trueProbDeath = overSampDat$eaDat$trueProbDeathNoNug
  SRSDat$eaDat$trueProbDeathNoNug = expit(logit(SRSDat$eaDat$trueProbDeathNoNug) - gamma*SRSDat$eaDat$urban)
  SRSDat$eaDat$trueProbDeath = SRSDat$eaDat$trueProbDeathNoNug
  overSampDat$eaDat$died = rbinom(nrow(overSampDat$eaDat), overSampDat$eaDat$numChildren, overSampDat$eaDat$trueProbDeathNoNug)
  SRSDat$eaDat$died = overSampDat$eaDat$died
  
  overSampDatTest$eaDat$trueProbDeathNoNug = expit(logit(overSampDatTest$eaDat$trueProbDeathNoNug) - gamma*overSampDatTest$eaDat$urban)
  overSampDatTest$eaDat$trueProbDeath = overSampDatTest$eaDat$trueProbDeathNoNug
  SRSDatTest$eaDat$trueProbDeathNoNug = expit(logit(SRSDatTest$eaDat$trueProbDeathNoNug) - gamma*SRSDatTest$eaDat$urban)
  SRSDatTest$eaDat$trueProbDeath = SRSDatTest$eaDat$trueProbDeathNoNug
  overSampDatTest$eaDat$died = rbinom(nrow(overSampDatTest$eaDat), overSampDatTest$eaDat$numChildren, overSampDatTest$eaDat$trueProbDeathNoNug)
  SRSDatTest$eaDat$died = overSampDatTest$eaDat$died
  for(i in 1:nsimBig) {
    overSampDat$clustDat[[i]]$trueProbDeathNoNug = expit(logit(overSampDat$clustDat[[i]]$trueProbDeathNoNug) - gamma*overSampDat$clustDat[[i]]$urban)
    overSampDat$clustDat[[i]]$trueProbDeath = overSampDat$clustDat[[i]]$trueProbDeathNoNug
    SRSDat$clustDat[[i]]$trueProbDeathNoNug = expit(logit(SRSDat$clustDat[[i]]$trueProbDeathNoNug) - gamma*SRSDat$clustDat[[i]]$urban)
    SRSDat$clustDat[[i]]$trueProbDeath = SRSDat$clustDat[[i]]$trueProbDeathNoNug
    overSampDatTest$clustDat[[i]]$trueProbDeathNoNug = expit(logit(overSampDatTest$clustDat[[i]]$trueProbDeathNoNug) - gamma*overSampDatTest$clustDat[[i]]$urban)
    overSampDatTest$clustDat[[i]]$trueProbDeath = overSampDatTest$clustDat[[i]]$trueProbDeathNoNug
    SRSDatTest$clustDat[[i]]$trueProbDeathNoNug = expit(logit(SRSDatTest$clustDat[[i]]$trueProbDeathNoNug) - gamma*SRSDatTest$clustDat[[i]]$urban)
    SRSDatTest$clustDat[[i]]$trueProbDeath = SRSDatTest$clustDat[[i]]$trueProbDeathNoNug
    
    overSampDat$clustDat[[i]]$died = overSampDat$eaDat$died[overSampClustDat$eaIs[,i]]
    SRSDat$clustDat[[i]]$died = SRSDat$eaDat$died[SRSClustDat$eaIs[,i]]
    overSampDatTest$clustDat[[i]]$died = overSampDatTest$eaDat$died[overSampClustDatTest$eaIs[,i]]
    SRSDatTest$clustDat[[i]]$died = SRSDatTest$eaDat$died[SRSClustDatTest$eaIs[,i]]
  }
  
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID02, "Big.RData"))
  overSampDat = overSampDatTest
  SRSDat = SRSDatTest
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID02, "TestBig.RData"))
  load(paste0("simDataMulti", dataID02, "Big.RData"))
  
  # Again, take only the first nsim simulations from the "big"" dataset
  overSampDat$clustDat = overSampDat$clustDat[1:nsim]
  SRSDat$clustDat = SRSDat$clustDat[1:nsim]
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID02, ".RData"))
  overSampDat = overSampDatTest
  SRSDat = SRSDatTest
  overSampDat$clustDat = overSampDat$clustDat[1:nsim]
  SRSDat$clustDat = SRSDat$clustDat[1:nsim]
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID02, "Test.RData"))
  
  # clustDat = SRSDat$clustDat[[1]]
  # clustDat = SRSDatTest$clustDat[[1]]
  # clustDat = overSampDat$clustDat[[1]]
  clustDat = overSampDatTest$clustDat[[1]]
  eaDat = overSampDat$eaDat
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = clustDat$numChildren
  obsCounts = clustDat$died
  
  pdf(paste0("figures/exampleOverSampTestSimulationNoUrban", dataID02, ".pdf"), width=8, height=8)
  par(mfrow =c(2, 2))
  zlim = c(0, quantile(c(eaDat$died/eaDat$numChildren, clustDat$died/clustDat$numChildren, 
                         eaDat$trueProbDeath), probs=.975))
  quilt.plot(eaDat$east, eaDat$north, eaDat$died/eaDat$numChildren, main="All Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, obsCounts/obsNs, main="Sample Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(eaDat$east, eaDat$north, eaDat$trueProbDeath, main="All True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, clustDat$trueProbDeath, main="Sample True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  dev.off()
  
  ##### do the same thing but now including only the intercept
  # reload the data 
  out = load(paste0("simDataMulti", dataID, "Big.RData"))
  
  # Now simulate the data without a cluster effect or urban effect but with the same underlying probability surface otherwise
  tausq = 0
  overSampDat$eaDat$trueProbDeathNoNug = expit(beta0)
  overSampDat$eaDat$trueProbDeath = overSampDat$eaDat$trueProbDeathNoNug
  SRSDat$eaDat$trueProbDeathNoNug = expit(beta0)
  SRSDat$eaDat$trueProbDeath = SRSDat$eaDat$trueProbDeathNoNug
  overSampDat$eaDat$died = rbinom(nrow(overSampDat$eaDat), overSampDat$eaDat$numChildren, overSampDat$eaDat$trueProbDeathNoNug)
  SRSDat$eaDat$died = overSampDat$eaDat$died
  
  overSampDatTest$eaDat$trueProbDeathNoNug = expit(beta0)
  overSampDatTest$eaDat$trueProbDeath = overSampDatTest$eaDat$trueProbDeathNoNug
  SRSDatTest$eaDat$trueProbDeathNoNug = expit(beta0)
  SRSDatTest$eaDat$trueProbDeath = SRSDatTest$eaDat$trueProbDeathNoNug
  overSampDatTest$eaDat$died = rbinom(nrow(overSampDatTest$eaDat), overSampDatTest$eaDat$numChildren, overSampDatTest$eaDat$trueProbDeathNoNug)
  SRSDatTest$eaDat$died = overSampDatTest$eaDat$died
  for(i in 1:nsimBig) {
    overSampDat$clustDat[[i]]$trueProbDeathNoNug = expit(beta0)
    overSampDat$clustDat[[i]]$trueProbDeath = overSampDat$clustDat[[i]]$trueProbDeathNoNug
    SRSDat$clustDat[[i]]$trueProbDeathNoNug = expit(beta0)
    SRSDat$clustDat[[i]]$trueProbDeath = SRSDat$clustDat[[i]]$trueProbDeathNoNug
    overSampDatTest$clustDat[[i]]$trueProbDeathNoNug = expit(beta0)
    overSampDatTest$clustDat[[i]]$trueProbDeath = overSampDatTest$clustDat[[i]]$trueProbDeathNoNug
    SRSDatTest$clustDat[[i]]$trueProbDeathNoNug = expit(beta0)
    SRSDatTest$clustDat[[i]]$trueProbDeath = SRSDatTest$clustDat[[i]]$trueProbDeathNoNug
    
    overSampDat$clustDat[[i]]$died = overSampDat$eaDat$died[overSampClustDat$eaIs[,i]]
    SRSDat$clustDat[[i]]$died = SRSDat$eaDat$died[SRSClustDat$eaIs[,i]]
    overSampDatTest$clustDat[[i]]$died = overSampDatTest$eaDat$died[overSampClustDatTest$eaIs[,i]]
    SRSDatTest$clustDat[[i]]$died = SRSDatTest$eaDat$died[SRSClustDatTest$eaIs[,i]]
  }
  
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID03, "Big.RData"))
  overSampDat = overSampDatTest
  SRSDat = SRSDatTest
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID03, "TestBig.RData"))
  load(paste0("simDataMulti", dataID03, "Big.RData"))
  
  # Again, take only the first nsim simulations from the "big"" dataset
  overSampDat$clustDat = overSampDat$clustDat[1:nsim]
  SRSDat$clustDat = SRSDat$clustDat[1:nsim]
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID03, ".RData"))
  overSampDat = overSampDatTest
  SRSDat = SRSDatTest
  overSampDat$clustDat = overSampDat$clustDat[1:nsim]
  SRSDat$clustDat = SRSDat$clustDat[1:nsim]
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID03, "Test.RData"))
  
  # clustDat = SRSDat$clustDat[[1]]
  # clustDat = SRSDatTest$clustDat[[1]]
  # clustDat = overSampDat$clustDat[[1]]
  clustDat = overSampDatTest$clustDat[[1]]
  eaDat = overSampDat$eaDat
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = clustDat$numChildren
  obsCounts = clustDat$died
  
  pdf(paste0("figures/exampleOverSampTestSimulationConstant", dataID03, ".pdf"), width=8, height=8)
  par(mfrow =c(2, 2))
  zlim = c(0, quantile(c(eaDat$died/eaDat$numChildren, clustDat$died/clustDat$numChildren, 
                         eaDat$trueProbDeath), probs=.975))
  quilt.plot(eaDat$east, eaDat$north, eaDat$died/eaDat$numChildren, main="All Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, obsCounts/obsNs, main="Sample Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(eaDat$east, eaDat$north, eaDat$trueProbDeath, main="All True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, clustDat$trueProbDeath, main="Sample True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  dev.off()
  
  setwd(wd)
  
  invisible(NULL)
}

# 


