source("setup.R")
index = commandArgs(trailingOnly = TRUE)
load("spdeCommandArgs.RData")
argList = spdeCommandArgs[[index]]
do.call("resultsSPDE", argList)