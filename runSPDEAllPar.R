source("setupParallel.R")
index = as.numeric(commandArgs(trailingOnly = TRUE))
load("spdeCommandArgs.RData")
argList = spdeCommandArgs[[index]]
do.call("resultsSPDEPar", argList)