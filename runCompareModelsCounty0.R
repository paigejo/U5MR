source("setup.R")
source("compareModels.R")
runCompareModels2(resultType="county", tausq=0, sampling="SRS", modelsI=1:7, printIEvery=10, nsim=10, saveResults=TRUE)
runCompareModels2(resultType="county", tausq=0, sampling="oversamp", modelsI=1:7, printIEvery=10, nsim=10, saveResults=TRUE)