# sets up global variables and libraries to run code in parallel
assign("doParallel", TRUE, envir=.GlobalEnv)

library(parallel)
library(doParallel)
library(foreach)
setwd("~/git/UM5R/")
source("setup.R")

# assign("cl", detectCores() - 1, envir=.GlobalEnv)
assign("cl", 7, envir=.GlobalEnv)
registerDoParallel(cl)
# stopCluster(cl)