# get all spde results
setwd("~/git/UM5R/")
source("setupParallel.R")
resultsSPDE(genCountLevel=TRUE)
resultsSPDE(genCountLevel=TRUE, tausq = .1^2)