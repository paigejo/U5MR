load("simDat.RData")
library(profvis)
source("setup.R")

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
nsim=100
set.seed(580252)
beta0 = -1.75
margVar = .15^2
tausq = .1^2
gamma = -1
# HHoldVar = .3^2
HHoldVar = 0
urbanOverSample = 2
effRange = 150

# there should be 1 true data set, but many simulated cluster samples
simulatedEAs = simDat2(kenyaEAs, clustDat=NULL, nsim=1, urbanOverSample=urbanOverSample,
                      beta0=beta0, margVar=margVar, tausq=tausq, gamma=gamma, HHoldVar=HHoldVar, 
                      effRange=effRange)

# simulate the cluster sampling and add to the data sets
overSampClustDat = simClusters3(kenyaEAs, urbanOverSample=urbanOverSample, nsim=nsim)
clustList = genAndreaFormatFromEAIs(simulatedEAs$eaDat, overSampClustDat$eaIs, overSampClustDat$sampleWeights)
overSampDat = list(eaDat=simulatedEAs$eaDat, clustDat=clustList)
SRSClustDat = simClusters3(kenyaEAs, urbanOverSample=1, nsim=nsim)
clustList = genAndreaFormatFromEAIs(simulatedEAs$eaDat, SRSClustDat$eaIs, SRSClustDat$sampleWeights)
SRSDat = list(eaDat=simulatedEAs$eaDat, clustDat=clustList) # the only thing different is the smapling of the clusters

# plot the first simulation of the over sampled and simple random sample data sets
clustDat = SRSDat$clustDat[[1]]
clustDat = overSampDat$clustDat[[1]]
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

save(overSampDat, SRSDat, file=paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                           round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, "urbanOver", 
                           round(urbanOverSample, 4), ".RData"))

# Now simulate the data without a cluster effect but with the same underlying probability surface otherwise
tausq = 0
overSampDat$eaDat$trueProbDeath = overSampDat$eaDat$trueProbDeathNoNug
SRSDat$eaDat$trueProbDeath = SRSDat$eaDat$trueProbDeathNoNug
overSampDat$eaDat$died = rbinom(nrow(overSampDat$eaDat), overSampDat$eaDat$numChildren, overSampDat$eaDat$trueProbDeathNoNug)
SRSDat$eaDat$died = rbinom(nrow(SRSDat$eaDat), SRSDat$eaDat$numChildren, SRSDat$eaDat$trueProbDeathNoNug)
for(i in 1:100) {
  overSampDat$clustDat[[i]]$trueProbDeath = overSampDat$clustDat[[i]]$trueProbDeathNoNug
  SRSDat$clustDat[[i]]$trueProbDeath = SRSDat$clustDat[[i]]$trueProbDeathNoNug
  
  overSampDat$clustDat[[i]]$died = overSampDat$eaDat$died[overSampClustDat$eaIs[,i]]
  SRSDat$clustDat[[i]]$died = SRSDat$eaDat$died[SRSClustDat$eaIs[,i]]
}
save(overSampDat, SRSDat, file=paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                                      round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, "urbanOver", 
                                      round(urbanOverSample, 4), ".RData"))

clustDat = SRSDat$clustDat[[1]]
clustDat = overSampDat$clustDat[[1]]
eaDat = overSampDat$eaDat
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