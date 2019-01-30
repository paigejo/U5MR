source("setup.R")
source("compareModels.R")
# runCompareModels2(test=TRUE, tausq=0, resultType="county", sampling="SRS", modelsI=1:7, printIEvery=1, nsim=10, saveResults=TRUE)
runCompareModels2(test=TRUE, tausq=0, resultType="county", sampling="oversamp", modelsI=1:7, printIEvery=1, nsim=10, saveResults=TRUE)