source("setup.R")
index = commandArgs(trailingOnly = TRUE)
load("directCommandArgs.RData")
print("index")
print(index)
argList = directCommandArgs[[index]]
do.call("getDirectNaive", argList)