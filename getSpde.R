# this script can be called as a batch job on the cluster, and should be used for generating spde results
setwd("~/git/UM5R/")
source("setup.R")
noNugget = resultsSPDE(genCountLevel=TRUE, tausq = 0)
nugget = resultsSPDE(genCountLevel=TRUE, tausq = .1^2)