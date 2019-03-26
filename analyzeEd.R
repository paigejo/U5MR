source("setup.R")
source("neonatalSimStudyWeighted.R")

# script for analyzing the women secondary education completion dataset

# first name elements of ed to be the same as the corresponding elements of the simulated datasets

##### first generate results for direct and naïve models
startIEd = 1
for(i in 1:nrow(ed)) {
  if(i %% 100 == 1)
    print(paste0("i: ", i))
  
  tmpDat = extendDataEd(ed[i,], v001 = i)
  
  # initialize the data frames on the first iteration only
  if(i == 1) {
    resEd = as.data.frame(matrix(nrow=sum(ed$n), ncol=ncol(tmpDat) + 1))
  }
  
  # append to the data frames
  names(resEd) = c(names(tmpDat), "regionRural")
  endIEd = startIEd + nrow(tmpDat) - 1
  resEd[startIEd:endIEd, 1:ncol(tmpDat)] = tmpDat
  
  # update row index
  startIEd = endIEd + 1
}

# add in RegionRural interaction
resEd$regionRural <- with(resEd, interaction(admin1, urban), drop=TRUE)

# save the resulting data frame
save(resEd, file=paste0("data4directEd.RData"))

# start analysis using naive approach and computing direct estimates
directEstEd = naiveEd = list()

# analyse the unstratified sampling scenario
education_obj2 = resEd
res2 = defineSurveyEd(education_obj2, 
                    stratVar=education_obj2$regionRural,
                    useSamplingWeights = TRUE)
directEstEd = res2
defineSurvey
resnA2 = run_naiveEd(education_obj2)

naiveEd = resnA2

# save results
save(directEstEd, naiveEd, file="resultsDirectNaiveEd.RData")

##### run Mercer et al. model
source("mercer.R")
tmpEd = mercer_u1m(directEstEd$logit.est, directEstEd$var.est, 
                    graph.path = "Kenyaadm1.graph")

resEd = data.frame(admin1=directEstEd$admin1,
                    est.mercer=expit(tmpEd$summary.linear.predictor$mean),
                    lower.mercer=tmpEd$summary.linear.predictor$"0.1quant",
                    upper.mercer=tmpEd$summary.linear.predictor$"0.9quant",
                    logit.est.mercer=tmpEd$summary.linear.predictor$mean, 
                    var.est.mercer=(tmpEd$summary.linear.predictor$sd)^2)
mercerEd = resEd

save(mercerEd, file=paste0("resultsMercerEd.RData"))

##### run BYM models
source("designBased.R")
runBYM2Ed(ed, includeUrbanRural = FALSE, includeCluster = FALSE)
runBYM2Ed(ed, includeUrbanRural = FALSE, includeCluster = TRUE)
runBYM2Ed(ed, includeUrbanRural = TRUE, includeCluster = FALSE)
runBYM2Ed(ed, includeUrbanRural = TRUE, includeCluster = TRUE)

##### run SPDE 
argList = list(list(clustDat = ed, includeClustEffect = FALSE, urbanEffect = FALSE), 
               list(clustDat = ed, includeClustEffect = FALSE, urbanEffect = TRUE), 
               list(clustDat = ed, includeClustEffect = TRUE, urbanEffect = FALSE), 
               list(clustDat = ed, includeClustEffect = TRUE, urbanEffect = TRUE))

for(i in 1:length(argList)) {
  args = argList[[i]]
  spdeResults = do.call("resultsSPDEed", args)
  includeClustEffect = args$includeClustEffect
  urbanEffect = args$urbanEffect
  fileName = paste0("resultsSPDEed_includeClustEffect", includeClustEffect, 
                    "_urbanEffect", urbanEffect, ".RData")
  save(spdeResults, file=fileName)
}







