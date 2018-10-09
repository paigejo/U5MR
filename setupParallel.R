# sets up global variables and libraries to run code in parallel
assign("doParallel", TRUE, envir=.GlobalEnv)

setwd("~/git/UM5R/")
library(parallel)
library(doParallel)
library(foreach)
source("setup.R")

# set up parallel environment, and generate a log
assign("cores", detectCores(), envir=.GlobalEnv)
assign("cl", makeCluster(cores[1]-1), envir=.GlobalEnv)
registerDoParallel(cl)
clusterEvalQ(cl, {setwd("~/git/UM5R/"); source("setup.R")})
writeLines(c(""), "log.txt")

# remember to run:
# stopCluster(cl)
# when you're done