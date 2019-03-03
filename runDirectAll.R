source("setup.R")
index = commandArgs(trailingOnly = TRUE)
load("directCommandArgs.RData")
argList = directCommandArgs[[index]]
do.call("getDirectNaive", argList)