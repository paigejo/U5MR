source("setup.R")
source("designBased.R")
runBYM(test=TRUE)
runBYM(test=TRUE, includeUrbanRural=FALSE)
runBYM(test=TRUE, includeCluster=FALSE)
runBYM(test=TRUE, includeUrbanRural=FALSE, includeCluster=FALSE)