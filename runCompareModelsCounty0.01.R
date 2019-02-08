source("setup.R")
source("compareModels.R")
runCompareModels2(resultType="county", sampling="SRS", modelsI=1:9, printIEvery=1, nsim=10, saveResults=TRUE)
runCompareModels2(resultType="county", sampling="oversamp", modelsI=1:9, printIEvery=1, nsim=10, saveResults=TRUE)