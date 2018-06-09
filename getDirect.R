source("neonatalSimStudyWeighted.R")

# load("simDataMulti.RData")
# load a different 1 of these depending on whether a cluster effect should be included 
# in the simulation of the data or not (tausq is the cluster effect variance)
load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOver2.RData")
# load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOver2.RData")

data4directSRS = list()
data4directoverSamp = list()

# extend the original simulated datasets to binary format
# where each child gets one row in the dataset.
for(j in 1:length(SRSDat[[2]])){
  resSRS = data.frame()
  resoverSamp = data.frame()
  for(i in 1:nrow(SRSDat[[2]][[j]])){
    resSRS = rbind(resSRS, extendData(SRSDat[[2]][[j]][i,], v001 = i))
    resoverSamp = rbind(resoverSamp, extendData(overSampDat[[2]][[j]][i,], v001 = i))
  }
  data4directSRS[[j]] = resSRS
  data4directoverSamp[[j]] = resoverSamp
}

#save(data4directSRS, data4directoverSamp, file="data4direct.RData")

# start analysis using naive approach and computing direct estimates
directEstSRS = naiveSRS = list()
directEstoverSamp = naiveoverSamp = list()

for(i in 1:100){
  print(i)
  # analyse the oversampled scenario
  childBirths_obj = data4directoverSamp[[i]]
  res = defineSurvey(childBirths_obj, 
                     stratVar=childBirths_obj$regionRural,
                     useSamplingWeights = TRUE)
  directEstoverSamp[[i]] = res
  
  resnA = run_naive(childBirths_obj)

  naiveoverSamp[[i]] = resnA
  
  # analyse the unstratified sampling scenario
  childBirths_obj2 = data4directSRS[[i]]
  res2 = defineSurvey(childBirths_obj2, 
                     stratVar=childBirths_obj2$regionRural,
                     useSamplingWeights = TRUE)
  directEstSRS[[i]] = res2
  
  resnA2 = run_naive(childBirths_obj2)
  
  naiveSRS[[i]] = resnA2
}

# save(directEstSRS, directEstoverSamp, naiveSRS, naiveoverSamp,file="resultsDirectNaiveTausq0.01.RData")
save(directEstSRS, directEstoverSamp, naiveSRS, naiveoverSamp,file="resultsDirectNaiveTausq0.RData")

