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

##### 
