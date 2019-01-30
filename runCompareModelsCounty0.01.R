source("setup.R")
source("compareModels.R")
runCompareModels2(resultType="county", sampling="SRS", modelsI=1:7, printIEvery=10, nsim=1, saveResults=TRUE)
runCompareModels2(resultType="county", sampling="oversamp", modelsI=1:7, printIEvery=10, nsim=1, saveResults=TRUE)