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
