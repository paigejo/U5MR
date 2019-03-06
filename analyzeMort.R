source("setup.R")

# script for analyzing the neonatal mortality dataset

##### first generate results for direct and naïve models
for(i in 1:nrow(mort)) {
  if(i %% 100 == 1)
    print(paste0("i: ", i))
  
  tmpDat = extendData(mort[i,], v001 = i)
  
  # initialize the data frames on the first iteration only
  if(i == 1) {
    resMort = as.data.frame(matrix(nrow=sum(MortDat[[2]][[j]]$numChildren), ncol=ncol(tmpDat) + 1))
  }
  
  # append to the data frames
  names(resMort) = c(names(tmpDat), "regionRural")
  endIMort = startIMort + nrow(tmpDat) - 1
  resMort[startIMort:endIMort, 1:ncol(tmpDat)] = tmpDat
  
  # update row index
  startIMort = endIMort + 1
}

# add in RegionRural interaction
resMort$regionRural <- with(resMort, interaction(admin1, urbanRural), drop=TRUE)

# save the resulting data frame
save(resMort, file=paste0("data4directMort.RData"))

# start analysis using naive approach and computing direct estimates
directEstMort = naiveMort = list()

# analyse the unstratified sampling scenario
childBirths_obj2 = resMort
res2 = defineSurvey(childBirths_obj2, 
                    stratVar=childBirths_obj2$regionRural,
                    useSamplingWeights = TRUE)
directEstMort = res2

resnA2 = run_naive(childBirths_obj2)

naiveMort = resnA2

# save results
save(directEstMort, naiveMort, file="resultsDirectNaiveMort.RData")

##### run Mercer et al. model
source("mercer.R")
tmpMort = mercer_u1m(directEstMort$logit.est, directEstMort$var.est, 
                    graph.path = "Kenyaadm1.graph")

resMort = data.frame(admin1=directEstMort$admin1,
                    u1m.mercer=expit(tmpMort$summary.linear.predictor$mean),
                    lower.mercer=tmpMort$summary.linear.predictor$"0.1quant",
                    upper.mercer=tmpMort$summary.linear.predictor$"0.9quant",
                    logit.est.mercer=tmpMort$summary.linear.predictor$mean, 
                    var.est.mercer=(tmpMort$summary.linear.predictor$sd)^2)
mercerMort = resMort

save(mercerMort, file=paste0("resultsMercerMort.RData"))

##### run BYM models
runBYM2Mort(mort, includeUrbanRural = FALSE, includeCluster = FALSE)
runBYM2Mort(mort, includeUrbanRural = FALSE, includeCluster = TRUE)
runBYM2Mort(mort, includeUrbanRural = FALSE, includeCluster = FALSE)
runBYM2Mort(mort, includeUrbanRural = TRUE, includeCluster = TRUE)

##### run SPDE 
# get prediction locations from population grid
# popGrid = makeInterpPopGrid()
load("popGrid.RData")
predCoords = cbind(popGrid$east, popGrid$north)
predUrban = popGrid$urban
# if(genEALevel) {
#   # Must predict at enumeration areas as well. Include enumeration areas as 
#   # first rows of prediction coordinates and prediction urban/rural
#   predCoords = rbind(cbind(eaDat$east, eaDat$north), predCoords)
#   predUrban = c(eaDat$urban, predUrban)
# }
# we only care about the probability, not counts, so not used except for the purposes 
# of calling inla:
# predNs = rep(25, nrow(predCoords))
predNs = rep(1, nrow(predCoords))

# get observations from dataset
obsCoords = cbind(mort$east, mort$north)
obsNs = mort$numChildren
obsCounts = mort$died
obsUrban = mort$urban



# fit model, get all predictions for each areal level and each posterior sample
fit = fitSPDEModel(obsCoords, obsNs=obsNs, obsCounts, obsUrban, predCoords, predNs=predNs, 
                   predUrban, genCountyLevel=TRUE, popGrid=popGrid, nPostSamples=nPostSamples, 
                   verbose = verbose, clusterEffect=includeClustEffect, 
                   int.strategy=int.strategy, genRegionLevel=genRegionLevel, 
                   keepPixelPreds=keepPixelPreds, genEALevel=genEALevel, 
                   urbanEffect=urbanEffect, link=1, predictionType=predictionType, 
                   exactAggregation=exactAggregation, genCountLevel=genCountLevel, 
                   eaDat=eaDat, truthByCounty=truthByCounty, truthByRegion=truthByRegion, 
                   truthByPixel=truthByPixel)
print(paste0("Fit completed: iteration ", i, "/", nsim))
countyPredMat = fit$countyPredMat
regionPredMat = fit$regionPredMat
pixelPredMat = fit$pixelPredMat
eaPredMat = fit$eaPredMat
eaMarginals = fit$eaMarginals
# fitSPDEModel3 = function(obsCoords, obsNs=rep(25, nrow(obsCoords)), obsCounts, obsUrban, predCoords, 
#                          predNs = rep(1, nrow(predCoords)), predUrban, clusterIndices, prior=NULL, 
#                          mesh=NULL, int.strategy="auto", strategy="laplace", 
#                          genCountyLevel=FALSE, popGrid=NULL, nPostSamples=100, kmRes=5, 
#                          counties=sort(unique(eaDat$admin1)), verbose=TRUE, genRegionLevel=FALSE, 
#                          regions=sort(unique(eaDat$region)), keepPixelPreds=FALSE, genEALevel=FALSE, 
#                          eaIndices=1:nrow(kenyaEAs), urbanEffect=TRUE, link=1, 
#                          predictionType=c("median", "mean"), eaDat=NULL, nSamplePixel=10, 
#                          truthByPixel=NULL, truthByCounty=NULL, truthByRegion=NULL, 
#                          truthByEa=NULL, clusterEffect=FALSE, significance=.8)









