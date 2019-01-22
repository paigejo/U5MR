# load("simDat.RData")
library(profvis)
library(logitnorm)
source("setup.R")
setwd("~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/")

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

nsim=100
nsimBig = 250
set.seed(580252)
beta0 = -1.75
margVar = .15^2
tausq = .1^2
gamma = -1
# HHoldVar = .3^2
HHoldVar = 0
# urbanOverSample = 2
effRange = 150
# urbanOverSamplefrac = 0.25
urbanOverSamplefrac = 0

# there should be 1 true data set, but many simulated cluster samples
# simulatedEAs = simDat2(kenyaEAs, clustDat=NULL, nsim=1, urbanOverSample=urbanOverSample,
#                        beta0=beta0, margVar=margVar, tausq=tausq, gamma=gamma, HHoldVar=HHoldVar, 
#                        effRange=effRange)
load("empiricalDistributions.RData")
simulatedEAs = simDatEmpirical(empiricalDistributions, kenyaEAs, clustDat=NULL, nsim=1, 
                               beta0=beta0, margVar=margVar, urbanOverSamplefrac=urbanOverSamplefrac, 
                               tausq=tausq, gamma=gamma, HHoldVar=HHoldVar, effRange=effRange)
kenyaEAs = simulatedEAs$eaDat
kenyaEAs$eaIs = 1:nrow(kenyaEAs)
kenyaEAsLong = kenyaEAs[rep(1:nrow(kenyaEAs), kenyaEAs$nHH),]

set.seed(1234)
# simulate the cluster sampling and add to the data sets
# overSampClustDat = simClusters3(kenyaEAs, urbanOverSample=urbanOverSample, nsim=nsim)
# overSampClustDat = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsim, NULL, 0, 25, 
#                                         urbanOverSamplefrac)
# temp = simClustersEmpirical(kenyaEAs, kenyaEAsLong, 2, NULL, urbanOverSamplefrac, verbose=FALSE, SRS = TRUE)
# out = profvis({temp = simClustersEmpirical(kenyaEAs, kenyaEAsLong, 5, NULL, urbanOverSamplefrac, verbose=FALSE)})
overSampClustDat = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsimBig, NULL, urbanOverSamplefrac, verbose=FALSE)
overSampClustDat = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsimBig, NULL, urbanOverSamplefrac, verbose=FALSE)
clustList = genAndreaFormatFromEAIs(simulatedEAs$eaDat, overSampClustDat$eaIs, overSampClustDat$sampleWeights)
overSampDat = list(eaDat=simulatedEAs$eaDat, clustDat=clustList)

overSampClustDatTest = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsimBig, NULL, urbanOverSamplefrac, fixedPerStrata=TRUE, nPerStrata=3, verbose=FALSE)
clustListTest = genAndreaFormatFromEAIs(simulatedEAs$eaDat, overSampClustDatTest$eaIs, overSampClustDatTest$sampleWeights)
overSampDatTest = list(eaDat=simulatedEAs$eaDat, clustDat=clustListTest)

# SRSClustDat = simClusters3(kenyaEAs, urbanOverSample=1, nsim=nsim)
# SRSClustDat = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsim, NULL, 0, 25)
SRSClustDat = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsimBig, NULL, SRS=TRUE, verbose=FALSE)
clustList = genAndreaFormatFromEAIs(simulatedEAs$eaDat, SRSClustDat$eaIs, SRSClustDat$sampleWeights)
SRSDat = list(eaDat=simulatedEAs$eaDat, clustDat=clustList) # the only thing different is the sampling of the clusters

SRSClustDatTest = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsimBig, NULL, fixedPerStrata=TRUE, nPerStrata=3, SRS=TRUE, verbose=FALSE)
clustListTest = genAndreaFormatFromEAIs(simulatedEAs$eaDat, SRSClustDatTest$eaIs, SRSClustDatTest$sampleWeights)
SRSDatTest = list(eaDat=simulatedEAs$eaDat, clustDat=clustListTest) # the only thing different is the sampling of the clusters

# plot the first simulation of the over sampled and simple random sample data sets
clustDat = SRSDat$clustDat[[1]]
# clustDat = overSampDat$clustDat[[1]]
eaDat = overSampDat$eaDat
pdf("figures/exampleSRSSimulation.pdf", width=8, height=8)
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

# save(overSampDat, SRSDat, file=paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
#                                       round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, "urbanOver", 
#                                       round(urbanOverSample, 4), ".RData"))
# save(overSampDat, SRSDat, file=paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
#                                       round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, "nUrbanClustersOver", 
#                                       round(numClustersUrbanOversamp, 4), ".RData"))
save(overSampDat, SRSDat, file=paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                                      round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, "urbanOverSamplefrac", 
                                      round(urbanOverSamplefrac, 4), "Big.RData"))
overSampDat = overSampDatTest
SRSDat = SRSDatTest
save(overSampDat, SRSDat, file=paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                                      round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, "urbanOverSamplefrac", 
                                      round(urbanOverSamplefrac, 4), "TestBig.RData"))
out = load(paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
            round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, "urbanOverSamplefrac", 
            round(urbanOverSamplefrac, 4), "Big.RData"))

# Now take only the first nsim simulations from the "big"" dataset
overSampDat$clustDat = overSampDat$clustDat[1:nsim]
SRSDat$clustDat = SRSDat$clustDat[1:nsim]
save(overSampDat, SRSDat, file=paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                                      round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, "urbanOverSamplefrac", 
                                      round(urbanOverSamplefrac, 4), ".RData"))
overSampDat = overSampDatTest
SRSDat = SRSDatTest
overSampDat$clustDat = overSampDat$clustDat[1:nsim]
SRSDat$clustDat = SRSDat$clustDat[1:nsim]
save(overSampDat, SRSDat, file=paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                                      round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, "urbanOverSamplefrac", 
                                      round(urbanOverSamplefrac, 4), "Test.RData"))

# reload the data
out = load(paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                  round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, "urbanOverSamplefrac", 
                  round(urbanOverSamplefrac, 4), "Big.RData"))

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
# save(overSampDat, SRSDat, file=paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
#                                       round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, "urbanOver", 
#                                       round(urbanOverSample, 4), ".RData"))
# save(overSampDat, SRSDat, file=paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
#                                       round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, "nUrbanClustersOver", 
#                                       round(numClustersUrbanOversamp, 4), ".RData"))
save(overSampDat, SRSDat, file=paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                                      round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, "urbanOverSamplefrac", 
                                      round(urbanOverSamplefrac, 4), "Big.RData"))
overSampDat = overSampDatTest
SRSDat = SRSDatTest
save(overSampDat, SRSDat, file=paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                                      round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, "urbanOverSamplefrac", 
                                      round(urbanOverSamplefrac, 4), "TestBig.RData"))
load(paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
            round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, "urbanOverSamplefrac", 
            round(urbanOverSamplefrac, 4), "Big.RData"))

# Again, take only the first nsim simulations from the "big"" dataset
overSampDat$clustDat = overSampDat$clustDat[1:nsim]
SRSDat$clustDat = SRSDat$clustDat[1:nsim]
save(overSampDat, SRSDat, file=paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                                      round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, "urbanOverSamplefrac", 
                                      round(urbanOverSamplefrac, 4), ".RData"))
overSampDat = overSampDatTest
SRSDat = SRSDatTest
overSampDat$clustDat = overSampDat$clustDat[1:nsim]
SRSDat$clustDat = SRSDat$clustDat[1:nsim]
save(overSampDat, SRSDat, file=paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                                      round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, "urbanOverSamplefrac", 
                                      round(urbanOverSamplefrac, 4), "Test.RData"))

# clustDat = SRSDat$clustDat[[1]]
# clustDat = SRSDatTest$clustDat[[1]]
# clustDat = overSampDat$clustDat[[1]]
clustDat = overSampDatTest$clustDat[[1]]
eaDat = overSampDat$eaDat
obsCoords = cbind(clustDat$east, clustDat$north)
obsNs = clustDat$numChildren
obsCounts = clustDat$died
wd = getwd()
setwd("~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/")
pdf("figures/exampleOverSampTestSimulationNoNug.pdf", width=8, height=8)
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

# check to make sure there are at least three clusters in each strata
# getStrata = function(dat) {
#   strata = dat$urban * 47 + match(dat$admin1, easpc$County) - 1 # subtract 1 since Mombasa is not rural
#   strata[strata  >= 46] = strata[strata  >= 46]-1 # subtract 1 since Nairobi is not rural
#   strata
# }
# test = sapply(SRSDat$clustDat, getStrata)
# temp = apply(test, 2, table)
# min(temp)

setwd("~/git/U5MR/")
