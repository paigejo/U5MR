library(colorspace)

# makes plots as well as parameter estimate tables for the example applications in the manuscript
makeAllPlots = function(dat=ed, meanRange, meanRange2, meanTicks, meanTicks2, meanTickLabels, meanTickLabels2, 
                        meanRangeSPDE, meanTicksSPDE, meanTickLabelsSPDE, sdRange, sdRange2, 
                        sdTicks, sdTicks2, sdTicksSPDE, sdTickLabels, sdTickLabels2, sdTickLabelsSPDE, 
                        meanRangeND, meanTicksND, meanTickLabelsND, sdRangeND, sdTicksND, sdTickLabelsND, 
                        meanRangeBYM2, meanTicksBYM2, meanTickLabelsBYM2, sdTicksBYM2, sdTickLabelsBYM2, 
                        varName="SCR", plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                        sdCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
                        ncols=29, relativeCols=makeRedGreenDivergingColors(ncols), plotUrbanMap=FALSE) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  
  print("generating data visualizations...")
  
  # plot the actual data
  png(file=paste0("figures/", resultNameRoot, "/clustersUrban", plotNameRoot, ".png"), width=500, height=500)
  par(oma=c( 0,0,0,3), mar=c(5.1, 4.1, 4.1, 6))
  urban = dat$urban
  plot(dat$lon[!urban], dat$lat[!urban], pch=19, col="green", main=paste0("Urban vs. rural clusters"), xlim=kenyaLonRange, 
       ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.2)
  points(dat$lon[urban], dat$lat[urban], pch=19, col="blue", cex=.2)
  # world(add=TRUE)
  plotMapDat(adm1)
  dev.off()
  
  # plot a map of urbanicity if requested (this can take ~10 minutes)
  if(plotUrbanMap) {
    # inside this if statement since it takes around ten minutes to run
    makeUrbanMap(kmres=1, savePlot=TRUE)
  }
  
  png(file=paste0("figures/", resultNameRoot, "/empirical", plotNameRoot, ".png"), width=500, height=500)
  par(oma=c( 0,0,0,3), mar=c(5.1, 4.1, 4.1, 6))
  quilt.plot(dat$lon, dat$lat, dat$y / dat$n, nx=150, ny=150, ylim=kenyaLatRange, xlim=kenyaLonRange, 
             xlab="Longitude", ylab="Latitude", main=paste0("Empirical ", varName), col=meanCols)
  # world(add=TRUE)
  plotMapDat(adm1)
  dev.off()
  
  png(file=paste0("figures/", resultNameRoot, "/empirical", plotNameRoot, "Logit.png"), width=500, height=500)
  par(oma=c( 0,0,0,3), mar=c(5.1, 4.1, 4.1, 6))
  ticks = pretty(seq(0, max(dat$y / dat$n), l=10), n=10)
  ticks = logit(ticks[-c(1, 11)])
  varRange = expit(range(ticks))
  # par( oma=c( 0,0,0,5)) # save some room for the legend
  plot(cbind(dat$lon, dat$lat), type="n", main=paste0("Kenya empirical ", varName), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  quilt.plot(cbind(dat$lon, dat$lat), logit(dat$y / dat$n), col=meanCols, 
             nx=100, ny=100, add.legend=FALSE, add=TRUE, zlim=logit(varRange))
  plotMapDat(adm1, lwd=.5)
  world(add=TRUE)
  # par( oma=c(0,0,0,2))
  image.plot(zlim=logit(varRange), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=meanCols, add = TRUE, axis.args=list(at=ticks, labels=expit(ticks)))
  dev.off()
  
  png(file=paste0("figures/", resultNameRoot, "/empirical", plotNameRoot, "Discrete.png"), width=500, height=500)
  par(oma=c( 0,0,0,3), mar=c(5.1, 4.1, 4.1, 6))
  totals = aggregate(dat$n, list(dat$admin1), FUN=sum)
  counts = aggregate(dat$y, list(dat$admin1), FUN=sum)
  plotMapDat(adm1, plotVar=counts$x / totals$x, new = TRUE, main=paste0("Empirical ", varName), cols=meanCols)
  dev.off()
  
  png(file=paste0("figures/", resultNameRoot, "/denominatorsDiscrete", plotNameRoot, ".png"), width=500, height=500)
  par(oma=c( 0,0,0,3), mar=c(5.1, 4.1, 4.1, 6))
  plotMapDat(adm1, plotVar=totals$x, new = TRUE, main=paste0("Sample size"), cols=meanCols)
  dev.off()
  
  png(file=paste0("figures/", resultNameRoot, "/empirical", varName, "DiscreteLogit.png"), width=500, height=500)
  par(oma=c( 0,0,0,3), mar=c(5.1, 4.1, 4.1, 6))
  totals = aggregate(dat$n, list(dat$admin1), FUN=sum)
  counts = aggregate(dat$y, list(dat$admin1), FUN=sum)
  plotMapDat(adm1, plotVar=counts$x / totals$x, new = TRUE, main=paste0("Empirical ", varName), col=meanCols, 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
  dev.off()
  
  png(file="figures/populationDiscrete.png", width=500, height=500)
  countyPops = poppc$popTotal
  sortI = sort(poppc$County, index.return=TRUE)$ix
  countyPops = countyPops[sortI]
  par(oma=c( 0,0,0,3), mar=c(5.1, 4.1, 4.1, 6))
  zlim = range(countyPops)
  plotMapDat(adm1, plotVar=countyPops, new = TRUE, main=paste0("Total population"), col=popCols, zlim=log(zlim), scaleFun=log, scaleFunInverse=exp)
  dev.off()
  
  # entirely for testing the extendData function
  out = load(paste0("data4direct", resultNameRoot, ".RData"))
  png(file=paste0("figures/", resultNameRoot, "/extended", varName, "Discrete.png"), width=500, height=500)
  testY = aggregate(resDat$y, list(resDat$admin1), FUN=sum)
  testN = aggregate(rep(1, nrow(resDat)), list(resDat$admin1), FUN=sum)
  plotMapDat(adm1, plotVar=testY$x / testN$x, new = TRUE, main=paste0("Extended ", varName), cols=meanCols)
  dev.off()
  
  # TODO: change scale limits
  thisPop = popGrid$popOrig
  totalPop = 43*10^6 # from DHS 2014 survey final report page 2
  thisPop = totalPop * thisPop / sum(thisPop) / 5^2 # population density per km^2
  png(file=paste0("figures/populationDensity.png"), width=500, height=500)
  par(oma=c( 0,0,0,3), mar=c(5.1, 4.1, 4.1, 6))
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=TeX("Population density (people/$km^2$)"), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  popRange = log(c(1, max(thisPop)))
  popTicks = log(pretty(exp(popRange)))
  popTicks = c(log(1), log(10), log(100), log(1000), popTicks[-c(1, 2, 4, 6)])
  tickLabels = as.character(exp(popTicks))
  tickLabels[7] = ""
  quilt.plot(cbind(popGrid$lon, popGrid$lat, log(thisPop)), nx=150, ny=150, add.legend=FALSE, add=TRUE, 
             zlim=range(popTicks), col=popCols)
  plotMapDat(adm1, lwd=.5)
  image.plot(zlim=popRange, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=popCols, add = TRUE, axis.args=list(at=popTicks, labels=tickLabels), legend.mar = 0)
  dev.off()
  
  png(file=paste0("figures/populationDensityDiscrete.png"), width=500, height=500)
  test = aggregate(thisPop*5^2, list(popGrid$admin1), FUN=sum)
  zlim = range(test$x)
  ticks = c(100000, 500000, 1000000, 2000000, 3000000, 4000000)
  plotMapDat(adm1, plotVar=test$x, new = TRUE, main=paste0("Integrated population density"), 
             col=popCols, zlim=log(zlim), scaleFun=log, scaleFunInverse=exp, ticks = ticks)
  dev.off()
  
  # plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE logit predictive SDs", typeText), ylim=kenyaLatRange, 
  #      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  # quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$sds), 
  #            nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(log(sdRange2)))
  # plotMapDat(adm1, lwd=.5)
  # points(dat$lon, dat$lat, pch=".")
  # image.plot(zlim=range(log(sdRange2)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
  #            col=cols, add = TRUE, axis.args=list(at=log(sdTicks2), labels=sdTickLabels2), legend.mar = 0)
  
  tab = totals
  tab$x = counts$x / totals$x
  tab
  
  ##### Naive and direct estimates
  print("plotting direct and naive estimates...")
  
  # save(directEstEd, naiveEd, file="resultsDirectNaiveEd.RData")
  out = load(paste0("resultsDirectNaive", resultNameRoot, ".RData"))
  plotName = paste0(resultNameRoot, "Preds")
  png(file=paste0("figures/", resultNameRoot, "/naive", plotName, ".png"), width=1000, height=1200)
  par(mfrow=c(2,2))
  zlim = range(c(expit(naiveResults$upper),expit(naiveResults$lower) ))
  plotMapDat(adm1, plotVar=naiveResults$est, new = TRUE, main=paste0("Naive ", varName, " estimates"), zlim=meanRange, cols=meanCols)
  plotMapDat(adm1, plotVar=sqrt(naiveResults$var.est), new = TRUE, main=paste0("Naive logit predictive SDs"), zlim=sdRange, cols=sdCols)
  plotMapDat(adm1, plotVar=expit(naiveResults$upper), new = TRUE, main=paste0("Naive 10th percentile"), zlim=meanRange, cols=meanCols)
  plotMapDat(adm1, plotVar=expit(naiveResults$lower), new = TRUE, main=paste0("Naive 90th percentile"), zlim=meanRange, cols=meanCols)
  dev.off()
  png(file=paste0("figures/", resultNameRoot, "/naive", plotName, "Logit.png"), width=1000, height=1200)
  par(mfrow=c(2,2), oma=c( 0,0,0,4))
  zlim = range(c(expit(naiveResults$upper),expit(naiveResults$lower) ))
  plotMapDat(adm1, plotVar=naiveResults$est, new = TRUE, main=paste0("Naive ", varName, " estimates"), 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
  plotMapDat(adm1, plotVar=sqrt(naiveResults$var.est), new = TRUE, main=paste0("Naive logit predictive SDs"), 
             zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, cols=sdCols)
  plotMapDat(adm1, plotVar=expit(naiveResults$upper), new = TRUE, main=paste0("Naive 10th percentile"), 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
  plotMapDat(adm1, plotVar=expit(naiveResults$lower), new = TRUE, main=paste0("Naive 90th percentile"), 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
  dev.off()
  png(file=paste0("figures/", resultNameRoot, "/direct", plotName, ".png"), width=800, height=1000)
  par(mfrow=c(2,2))
  zlim = range(c(expit(directEstResults$upper),expit(directEstResults$lower) ))
  plotMapDat(adm1, plotVar=directEstResults$est, new = TRUE, main=paste0("Direct ", varName, " estimates"), 
             zlim=meanRange, cols=meanCols)
  plotMapDat(adm1, plotVar=sqrt(directEstResults$var.est), new = TRUE, main=paste0("Direct logit predictive SDs"), 
             zlim=sdRange, cols=sdCols)
  plotMapDat(adm1, plotVar=expit(directEstResults$upper), new = TRUE, main=paste0("Direct 10th percentile"), 
             zlim=meanRange, cols=meanCols)
  plotMapDat(adm1, plotVar=expit(directEstResults$lower), new = TRUE, main=paste0("Direct 90th percentile"), 
             zlim=meanRange, cols=meanCols)
  dev.off()
  png(file=paste0("figures/", resultNameRoot, "/direct", plotName, "Logit.png"), width=800, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,0,4))
  plotMapDat(adm1, plotVar=directEstResults$est, new = TRUE, main=paste0("Direct ", varName, " estimates"), 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
  plotMapDat(adm1, plotVar=sqrt(directEstResults$var.est), new = TRUE, main=paste0("Direct logit predictive SDs"), 
             zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, cols=sdCols)
  plotMapDat(adm1, plotVar=expit(directEstResults$upper), new = TRUE, main=paste0("Direct 10th percentile"), 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
  plotMapDat(adm1, plotVar=expit(directEstResults$lower), new = TRUE, main=paste0("Direct 90th percentile"), 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
  dev.off()
  
  ##### Mercer et al. estimates
  print("plotting mercer estimates...")
  
  # save(mercerEd, file=paste0("resultsMercerEd.RData"))
  plotName = paste0("mercer", resultNameRoot, "Preds")
  out = load(paste0("resultsMercer", resultNameRoot, ".RData" ))
  zlim = range(c(expit(mercerResults$lower.mercer),expit(mercerResults$upper.mercer)))
  png(file=paste0("figures/", resultNameRoot, "/", plotName, ".png"), width=800, height=1000)
  par(mfrow=c(2,2))
  plotMapDat(adm1, plotVar=mercerResults$est.mercer, new = TRUE, main=paste0("Mercer et al. ", varName, " estimates"), 
             zlim=meanRange, cols=meanCols)
  plotMapDat(adm1, plotVar=sqrt(mercerResults$var.est.mercer), new = TRUE, main=paste0("Mercer et al. logit predictive SDs"), 
             zlim=sdRange, cols=sdCols)
  plotMapDat(adm1, plotVar=expit(mercerResults$lower.mercer), new = TRUE, main=paste0("Mercer et al. 10th percentile"), 
             zlim=meanRange, cols=meanCols)
  plotMapDat(adm1, plotVar=expit(mercerResults$upper.mercer), new = TRUE, main=paste0("Mercer et al. 90th percentile"), 
             zlim=meanRange, cols=meanCols)
  dev.off()
  
  png(file=paste0("figures/", resultNameRoot, "/", plotName, "Logit.png"), width=800, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,0,4))
  plotMapDat(adm1, plotVar=mercerResults$est.mercer, new = TRUE, main=paste0("Mercer et al. ", varName, " estimates"), 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
  plotMapDat(adm1, plotVar=sqrt(mercerResults$var.est.mercer), new = TRUE, main=paste0("Mercer et al. logit predictive SDs"), 
             zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, cols=sdCols)
  plotMapDat(adm1, plotVar=expit(mercerResults$lower.mercer), new = TRUE, main=paste0("Mercer et al. 10th percentile"), 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
  plotMapDat(adm1, plotVar=expit(mercerResults$upper.mercer), new = TRUE, main=paste0("Mercer et al. 90th percentile"), 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
  dev.off()
  
  zlim = range(c(expit(mercerResults$lower.mercer),expit(mercerResults$upper.mercer)))
  png(file=paste0("figures/", resultNameRoot, "/", plotName, "2.png"), width=800, height=1000)
  par(mfrow=c(2,2))
  plotMapDat(adm1, plotVar=mercerResults$est.mercer, new = TRUE, main=paste0("Mercer et al. ", varName, " estimates"), 
             zlim=meanRange2, cols=meanCols)
  plotMapDat(adm1, plotVar=sqrt(mercerResults$var.est.mercer), new = TRUE, main=paste0("Mercer et al. logit predictive SDs"), 
             zlim=sdRange2, cols=sdCols)
  plotMapDat(adm1, plotVar=expit(mercerResults$lower.mercer), new = TRUE, main=paste0("Mercer et al. 10th percentile"), 
             zlim=meanRange2, cols=meanCols)
  plotMapDat(adm1, plotVar=expit(mercerResults$upper.mercer), new = TRUE, main=paste0("Mercer et al. 90th percentile"), 
             zlim=meanRange2, cols=meanCols)
  dev.off()
  
  png(file=paste0("figures/", resultNameRoot, "/", plotName, "Logit2.png"), width=800, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,0,4))
  plotMapDat(adm1, plotVar=mercerResults$est.mercer, new = TRUE, main=paste0("Mercer et al. ", varName, " estimates"), 
             zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, cols=meanCols)
  plotMapDat(adm1, plotVar=sqrt(mercerResults$var.est.mercer), new = TRUE, main=paste0("Mercer et al. logit predictive SDs"), 
             zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2, cols=sdCols)
  plotMapDat(adm1, plotVar=expit(mercerResults$lower.mercer), new = TRUE, main=paste0("Mercer et al. 10th percentile"), 
             zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, cols=meanCols)
  plotMapDat(adm1, plotVar=expit(mercerResults$upper.mercer), new = TRUE, main=paste0("Mercer et al. 90th percentile"), 
             zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, cols=meanCols)
  dev.off()
  
  ##### BYM2 estimates
  print("plotting BYM2 estimates...")
  
  # save(file = paste0('bym2EdUrbRur',includeUrbanRural, 'Cluster', includeCluster, '.RData'), 
  #      designRes = designRes)
  # if(includeCluster) {
  #   save(file = paste0('bym2EdUrbRur',includeUrbanRural, 'Cluster', includeCluster, 'debiased.RData'), 
  #        designRes = designRes)
  # }
  argList = list(list(includeUrbanRural = FALSE, includeCluster = FALSE), 
                 list(includeUrbanRural = FALSE, includeCluster = TRUE), 
                 list(includeUrbanRural = TRUE, includeCluster = FALSE), 
                 list(includeUrbanRural = TRUE, includeCluster = TRUE))
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    zlim = range(c(expit(designRes$predictions$Q10),expit(designRes$predictions$Q90)))
    png(file=paste0("figures/", resultNameRoot, "/", nameRoot, ".png"), width=800, height=1000)
    par(mfrow=c(2,2))
    plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", varName, " estimates", typeText), 
               zlim=meanRange, cols=meanCols)
    plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 logit predictive SDs", typeText), 
               zlim=sdRange, cols=sdCols)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 10th percentile", typeText), 
               zlim=meanRange, cols=meanCols)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 90th percentile", typeText), 
               zlim=meanRange, cols=meanCols)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/", nameRoot, "Logit.png"), width=800, height=1000)
    par(mfrow=c(2,2), oma=c( 0,0,0,4))
    plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", varName, " estimates", typeText), 
               zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
    plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 logit predictive SDs", typeText), 
               zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, cols=sdCols)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 10th percentile", typeText), 
               zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 90th percentile", typeText), 
               zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/", nameRoot, "2.png"), width=800, height=1000)
    par(mfrow=c(2,2))
    plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", varName, " estimates", typeText), cols=meanCols, zlim=meanRange2)
    plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 logit predictive SDs", typeText), cols=sdCols, zlim=sdRange2)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 10th percentile", typeText), cols=meanCols, zlim=meanRange2)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 90th percentile", typeText), cols=meanCols, zlim=meanRange2)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/", nameRoot, "Logit2.png"), width=800, height=1000)
    par(mfrow=c(2,2), oma=c( 0,0,0,4))
    plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", varName, " estimates", typeText), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
    plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 logit predictive SDs", typeText), cols=sdCols, zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 10th percentile", typeText), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 90th percentile", typeText), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
    dev.off()
    
    if(includeCluster) {
      # also gather and plot the debiased results
      nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
      out = load(paste0(nameRoot, '.RData'))
      
      urbanText = ifelse(includeUrban, "", "noUrb")
      clusterText = ifelse(includeCluster, "", "NoClust")
      both = includeUrban && includeUrban
      if(both)
        debiasedText = "Debiased"
      else
        debiasedText = "debiased"
      typeText = paste0(" ", urbanText, clusterText, debiasedText)
      
      png(file=paste0("figures/", resultNameRoot, "/", nameRoot, ".png"), width=800, height=1000)
      par(mfrow=c(2,2))
      plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", varName, " estimates", typeText), cols=meanCols, zlim=meanRange)
      plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 logit predictive SDs", typeText), cols=sdCols, zlim=sdRange)
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 10th percentile", typeText), cols=meanCols, zlim=meanRange)
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 90th percentile", typeText), cols=meanCols, zlim=meanRange)
      dev.off()
      
      png(file=paste0("figures/", resultNameRoot, "/", nameRoot, "Logit.png"), width=800, height=1000)
      par(mfrow=c(2,2), oma=c( 0,0,0,4))
      plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", varName, " estimates", typeText), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
      plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 logit predictive SDs", typeText), cols=sdCols, zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels)
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 10th percentile", typeText), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 90th percentile", typeText), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
      dev.off()
      
      png(file=paste0("figures/", resultNameRoot, "/", nameRoot, "2.png"), width=800, height=1000)
      par(mfrow=c(2,2))
      plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", varName, " estimates", typeText), cols=meanCols, zlim=meanRange2)
      plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 logit predictive SDs", typeText), cols=sdCols, zlim=sdRange2)
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 10th percentile", typeText), cols=meanCols, zlim=meanRange2)
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 90th percentile", typeText), cols=meanCols, zlim=meanRange2)
      dev.off()
      
      png(file=paste0("figures/", resultNameRoot, "/", nameRoot, "Logit2.png"), width=800, height=1000)
      par(mfrow=c(2,2), oma=c( 0,0,0,4))
      plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", varName, " estimates", typeText), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
      plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 logit predictive SDs", typeText), cols=sdCols, zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2)
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 10th percentile", typeText), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 90th percentile", typeText), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
      dev.off()
    }
  }
  
  ##### SPDE estimates
  print("plotting SPDE estimates...")
  
  argList = list(list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    zlim = range(c(spdeResults$resultsCounty$lower,spdeResults$resultsCounty$upper))
    
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, ".png"), width=800, height=1000)
    par(mfrow=c(2,2))
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE ", varName, " estimates", typeText), cols=meanCols, zlim=meanRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE logit predictive SDs", typeText), cols=sdCols, zlim=sdRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE 10th percentile", typeText), cols=meanCols, zlim=meanRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE 90th percentile", typeText), cols=meanCols, zlim=meanRange)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "Logit.png"), width=800, height=1000)
    par(mfrow=c(2,2), oma=c( 0,0,0,4))
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE ", varName, " estimates", typeText), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE logit predictive SDs", typeText), cols=sdCols, zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE 10th percentile", typeText), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE 90th percentile", typeText), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "2.png"), width=800, height=1000)
    par(mfrow=c(2,2))
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE ", varName, " estimates", typeText), cols=meanCols, zlim=meanRange2)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE logit predictive SDs", typeText), cols=sdCols, zlim=sdRange2)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE 10th percentile", typeText), cols=meanCols, zlim=meanRange2)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE 90th percentile", typeText), cols=meanCols, zlim=meanRange2)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "Logit2.png"), width=800, height=1000)
    par(mfrow=c(2,2), oma=c( 0,0,0,4))
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE ", varName, " estimates", typeText), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE logit predictive SDs", typeText), cols=sdCols, zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE 10th percentile", typeText), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE 90th percentile", typeText), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
    dev.off()
    
    ## plot continuous prediction surface
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, ".png"), width=800, height=1000)
    par(mfrow=c(2,2))
    quilt.plot(popGrid$lon, popGrid$lat, spdeResults$resultsPixel$pred, nx=200, ny=200, col=meanCols, zlim=meanRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE ", varName, " estimates", typeText), cols=meanCols, zlim=meanRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE logit predictive SDs", typeText), cols=sdCols, zlim=sdRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE 10th percentile", typeText), cols=meanCols, zlim=meanRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE 90th percentile", typeText), cols=meanCols, zlim=meanRange)
    dev.off()
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "Logit.png"), width=800, height=1000)
    par(mfrow=c(2,2), oma=c( 0,0,0,4))
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE ", varName, " estimates", typeText), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE logit predictive SDs", typeText), cols=sdCols, zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE 10th percentile", typeText), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE 90th percentile", typeText), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "ContinuousLogit.png"), width=800, height=1000)
    par(mfrow=c(2,2), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", varName, " estimates", typeText), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$pred), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRange)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE logit predictive SDs", typeText), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
    quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$sds), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=sdCols, zlim=range(log(sdRange)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(log(sdRange)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=sdCols, add = TRUE, axis.args=list(at=log(sdTicks), labels=sdTickLabels), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE 10th percentile", typeText), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$lower), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRange)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE 90th percentile", typeText), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$upper), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRange)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels), legend.mar = 5)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "ContinuousLogit2.png"), width=800, height=1000)
    par(mfrow=c(2,2), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", varName, " estimates", typeText), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$pred), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange2)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRange2)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks2), labels=meanTickLabels2), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE logit predictive SDs", typeText), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
    quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$sds), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=sdCols, zlim=range(log(sdRange2)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(log(sdRange2)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=sdCols, add = TRUE, axis.args=list(at=log(sdTicks2), labels=sdTickLabels2), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE 10th percentile", typeText), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$lower), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange2)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRange2)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks2), labels=meanTickLabels2), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE 90th percentile", typeText), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$upper), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange2)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRange2)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks2), labels=meanTickLabels2), legend.mar = 5)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "ContinuousLogitSelf.png"), width=800, height=1000)
    par(mfrow=c(2,2), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", varName, " estimates", typeText), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$pred), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE logit predictive SDs", typeText), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
    quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$sds), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=sdCols, zlim=range(log(sdRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(log(sdRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=sdCols, add = TRUE, axis.args=list(at=log(sdTicksSPDE), labels=sdTickLabelsSPDE), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE 10th percentile", typeText), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$lower), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE 90th percentile", typeText), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$upper), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 5)
    dev.off()
  }
  
  ##### now put the predictions from each together in the same plot when relevant
  # 1: put naive and direct together with standard deviations and credible intervals
  ## make two by four plot with these two models and lower, mean, and upper bounds, and SDs
  # 2: put values from the smooth models including all effects together:
  ## mercer, BYM2 with urban and cluster effects, SPDE with urban and cluster effects
  ## make three by four plot with these three models and lower, mean, and upper bounds, and SDs
  # 3: plot test spde models together
  ## 4 x 4 plot showing estimates, sds, 10% and 90% CIs
  # 4: plot relative differences between the models, relative to the BYM2 model including urban and cluster effects that has been debiased
  ## 2 x 2 plot
  # 5: plot the BYM2 models together
  ## 4 x 6 plot
  # 6: plot relative differences between the BYM2 models
  ## 3 x 2 plot (one of the subplots will be empty)
  # 7: plot relative differences between the SPDE models on a continuous scale
  ## 2 x 2 plot, all relative to the model with both urban and cluster effects
  # 8: plot relative differences between the SPDE models on a discrete scale
  ## 2 x 2 plot
  # 9: make pairs plot of the full models with respect to each other
  
  ## Plot 1: direct and naive models
  print("plotting direct and naive estimates together...")
  
  out = load(paste0("resultsDirectNaive", resultNameRoot, ".RData" ))
  
  png(file=paste0("figures/", resultNameRoot, "/fullDirectNaive", plotNameRoot, ".png"), width=800, height=1200)
  par(mfrow=c(4,2), oma=c( 0,0,0,4))
  plotMapDat(adm1, plotVar=naiveResults$est, new = TRUE, main=paste0("Naive ", varName, " estimates"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND)
  plotMapDat(adm1, plotVar=directEstResults$est, new = TRUE, main=paste0("Direct ", varName, " estimates"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND)
  plotMapDat(adm1, plotVar=sqrt(naiveResults$var.est), new = TRUE, main=paste0("Naive logit predictive SDs"), cols=sdCols, zlim=log(sdRangeND), scaleFun=log, scaleFunInverse=exp, ticks=sdTicksND, tickLabels=sdTickLabelsND)
  plotMapDat(adm1, plotVar=sqrt(directEstResults$var.est), new = TRUE, main=paste0("Direct logit predictive SDs"), cols=sdCols, zlim=log(sdRangeND), scaleFun=log, scaleFunInverse=exp, ticks=sdTicksND, tickLabels=sdTickLabelsND)
  plotMapDat(adm1, plotVar=expit(naiveResults$upper), new = TRUE, main=paste0("Naive 10th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND)
  plotMapDat(adm1, plotVar=expit(directEstResults$upper), new = TRUE, main=paste0("Direct 10th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND)
  plotMapDat(adm1, plotVar=expit(naiveResults$lower), new = TRUE, main=paste0("Naive 90th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND)
  plotMapDat(adm1, plotVar=expit(directEstResults$lower), new = TRUE, main=paste0("Direct 90th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND)
  dev.off()
  
  ## Plot 2: mercer, BYM2 with urban and cluster effects, SPDE with urban and cluster effects
  print("plotting smoothed estimates together...")
  
  png(file=paste0("figures/", resultNameRoot, "/fullSmoothed", plotNameRoot, ".png"), width=1000, height=1200)
  par(mfrow=c(4,3), oma=c( 0,0,0,4))
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
  out = load(paste0(nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeUrban
  debiasedText = "Debiased"
  typeTextBYM = paste0(" ", urbanText, clusterText, debiasedText)
  
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeUrban
  notBothText = ifelse(both, "", " ")
  typeTextSPDE = paste0(notBothText, urbanText, clusterText)
  
  plotMapDat(adm1, plotVar=mercerResults$est.mercer, new = TRUE, main=paste0("Mercer et al. ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", varName, " estimates", typeTextBYM), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE ", varName, " estimates", typeTextSPDE), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
  
  plotMapDat(adm1, plotVar=sqrt(mercerResults$var.est.mercer), new = TRUE, main=paste0("Mercer et al. logit predictive SDs"), cols=sdCols, zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2)
  plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 logit predictive SDs", typeTextBYM), cols=sdCols, zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE logit predictive SDs", typeTextSPDE), cols=sdCols, zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2)
  
  plotMapDat(adm1, plotVar=expit(mercerResults$lower.mercer), new = TRUE, main=paste0("Mercer et al. 10th percentile"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 10th percentile", typeTextBYM), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE 10th percentile", typeTextSPDE), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
  
  plotMapDat(adm1, plotVar=expit(mercerResults$upper.mercer), new = TRUE, main=paste0("Mercer et al. 90th percentile"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 90th percentile", typeTextBYM), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE 90th percentile", typeTextSPDE), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
  dev.off()
  
  ## Plot 3: plot all of the SPDE plots together (4 x 4 plot)
  print("plotting SPDE models together...")
  
  argList = list(list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "ContinuousTogetherSelf.png"), width=1400, height=1400)
  par(mfrow=c(4,4), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", varName, " estimates", typeText), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$pred), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 0)
  }
  
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE logit predictive SDs", typeText), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
    quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$sds), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=sdCols, zlim=range(log(sdRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(log(sdRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=sdCols, add = TRUE, axis.args=list(at=log(sdTicksSPDE), labels=sdTickLabelsSPDE), legend.mar = 0)
  }
  
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE 10th percentile", typeText), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$lower), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 0)
  }
  
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE 90th percentile", typeText), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$upper), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 5)
  }
  dev.off()
  
  ##### Plot 4: plot relative differences between the models
  print("plotting relative differences between models...")
  
  ## all together relative to the bym2
  # estimates
  # library("colorspace")
  # pal <-choose_palette()
  # cols=diverging_hcl(101, h1=265, h2=101, c1=100, l1=50, l2=92, p1=0.6, p2=1.5)
  # cols=
  
  nPerSide = (ncols - 1) / 2
  # cols=diverging_hcl(ncols, h1=265, h2=17, c1=100, l1=50, l2=92, p1=0.6)
  cols=relativeCols
  png(file=paste0("figures/", resultNameRoot, "/fullRelative", plotNameRoot, ".png"), width=800, height=1000)
  par(mfrow=c(2, 2), oma=c( 0,0,0,4))
  
  zlim = range(c((plotVar=naiveResults$est - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), 
                 (directEstResults$est - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), 
                 (mercerResults$est.mercer - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), 
                 (spdeResults$resultsCounty$pred - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean)))
  lengthUp = zlim[2] - 0
  lengthDown = 0 - zlim[1]
  if(lengthUp >= lengthDown) {
    # shorten the blue part of the color scale
    propDown = lengthDown / lengthUp
    propUp = 1
  } else {
    propDown = 1
    propUp = lengthUp / lengthDown
  }
  # cols = pal(101) # keep an extra color in the middle for 0
  numDown = round(nPerSide * propDown)
  numUp = round(nPerSide * propUp)
  cols = c(cols[(nPerSide + 1-numDown):nPerSide], cols[nPerSide + 1], cols[(nPerSide + 2):(nPerSide + 1 + numUp)])
  plotMapDat(adm1, plotVar=(plotVar=naiveResults$est - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), 
             new = TRUE, main=paste0("Naive relative to BYM2 ", varName, " estimates"), zlim=zlim, col=cols)
  plotMapDat(adm1, plotVar=(directEstResults$est - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), 
             new = TRUE, main=paste0("Direct relative to BYM2 ", varName, " estimates"), zlim=zlim, col=cols)
  plotMapDat(adm1, plotVar=(mercerResults$est.mercer - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), 
             new = TRUE, main=paste0("Mercer et al. relative to BYM2 ", varName, " estimates"), zlim=zlim, col=cols)
  plotMapDat(adm1, plotVar=(spdeResults$resultsCounty$pred - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), 
             new = TRUE, main=paste0("SPDE relative to BYM2 ", varName, " estimates", typeTextSPDE), zlim=zlim, col=cols)
  dev.off()
  
  ## Plot 5: plot the BYM2 models together (four by six plot)
  argList = list(list(includeUrbanRural = FALSE, includeCluster = FALSE), 
                 list(includeUrbanRural = FALSE, includeCluster = TRUE), 
                 list(includeUrbanRural = TRUE, includeCluster = FALSE), 
                 list(includeUrbanRural = TRUE, includeCluster = TRUE))
  png(file=paste0("figures/", resultNameRoot, "/predsBYM2togetherSelf", resultNameRoot, "Logit.png"), width=1500, height=1000)
  par(mfrow=c(4,6))
  
  # plot the estimates
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", varName, " estimates", typeText), 
               zlim=logit(meanRangeBYM2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksBYM2, tickLabels=meanTickLabelsBYM2, 
               cols=meanCols)
    
    if(includeCluster) {
      # also gather and plot the debiased results
      nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
      out = load(paste0(nameRoot, '.RData'))
      
      urbanText = ifelse(includeUrban, "", "noUrb")
      clusterText = ifelse(includeCluster, "", "NoClust")
      both = includeUrban && includeUrban
      if(both)
        debiasedText = "Debiased"
      else
        debiasedText = "debiased"
      typeText = paste0(" ", urbanText, clusterText, debiasedText)
      
      
      plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", varName, " estimates", typeText), 
                 cols=meanCols, zlim=logit(meanRangeBYM2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksBYM2, tickLabels=meanTickLabelsBYM2)
    }
  }
  
  # plot the BYM2 standard deviations together
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 logit predictive SDs", typeText), 
               zlim=log(sdRangeBYM2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicksBYM2, tickLabels=sdTickLabelsBYM2, cols=sdCols)
    
    if(includeCluster) {
      # also gather and plot the debiased results
      nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
      out = load(paste0(nameRoot, '.RData'))
      
      urbanText = ifelse(includeUrban, "", "noUrb")
      clusterText = ifelse(includeCluster, "", "NoClust")
      both = includeUrban && includeUrban
      if(both)
        debiasedText = "Debiased"
      else
        debiasedText = "debiased"
      typeText = paste0(" ", urbanText, clusterText, debiasedText)
      
      plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 logit predictive SDs", typeText), 
                 cols=sdCols, zlim=log(sdRangeBYM2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicksBYM2, 
                 tickLabels=sdTickLabelsBYM2)
    }
  }
  
  # plot the BYM2 tenth percentiles together
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 10th percentile", typeText), 
               zlim=logit(meanRangeBYM2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksBYM2, tickLabels=meanTickLabelsBYM2, cols=meanCols)
    
    if(includeCluster) {
      # also gather and plot the debiased results
      nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
      out = load(paste0(nameRoot, '.RData'))
      
      urbanText = ifelse(includeUrban, "", "noUrb")
      clusterText = ifelse(includeCluster, "", "NoClust")
      both = includeUrban && includeUrban
      if(both)
        debiasedText = "Debiased"
      else
        debiasedText = "debiased"
      typeText = paste0(" ", urbanText, clusterText, debiasedText)
      
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 10th percentile", typeText), 
                 cols=meanCols, zlim=logit(meanRangeBYM2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksBYM2, 
                 tickLabels=meanTickLabelsBYM2)
    }
  }
  
  # plot the BYM2 ninetieth percentiles together
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 90th percentile", typeText), 
               zlim=logit(meanRangeBYM2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksBYM2, tickLabels=meanTickLabelsBYM2, cols=meanCols)
    
    if(includeCluster) {
      # also gather and plot the debiased results
      nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
      out = load(paste0(nameRoot, '.RData'))
      
      urbanText = ifelse(includeUrban, "", "noUrb")
      clusterText = ifelse(includeCluster, "", "NoClust")
      both = includeUrban && includeUrban
      if(both)
        debiasedText = "Debiased"
      else
        debiasedText = "debiased"
      typeText = paste0(" ", urbanText, clusterText, debiasedText)
      
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 90th percentile", typeText), 
                 cols=meanCols, zlim=logit(meanRangeBYM2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksBYM2, 
                 tickLabels=meanTickLabelsBYM2)
    }
  }
  dev.off()
  
  ## Plot 6: plot the BYM2 model relative differences versus the full BYM2 model
  # 2 x 3 plot. Columns: debiased cluster, cluster, no cluster. Rows: urban, no urban
  argList = list(list(includeUrbanRural = TRUE, includeCluster = TRUE), 
                 list(includeUrbanRural = TRUE, includeCluster = FALSE), 
                 list(includeUrbanRural = FALSE, includeCluster = TRUE), 
                 list(includeUrbanRural = FALSE, includeCluster = FALSE))
  png(file=paste0("figures/", resultNameRoot, "/fullRelativeBYM2", resultNameRoot, ".png"), width=1200, height=800)
  par(mfrow=c(2,3))
  
  # before the plotting begins, determine the range of relative differences
  zlim=c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    
    if(includeCluster) {
      nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
      out = load(paste0(nameRoot, '.RData'))
      
      if(i == 1) {
        # get the values relative to which other predictions will be compared
        expected = expit(designRes$predictions$mean)
      }
      else {
        vals = (expit(designRes$predictions$mean) - expected) / expected
        zlim = range(c(zlim, vals))
      }
    }
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    vals = (expit(designRes$predictions$mean) - expected) / expected
    zlim = range(c(zlim, vals))
  }
  
  # now shift the color scales so that it center occurs at 0
  cols = relativeCols
  lengthUp = zlim[2] - 0
  lengthDown = 0 - zlim[1]
  if(lengthUp >= lengthDown) {
    # shorten the blue part of the color scale
    propDown = lengthDown / lengthUp
    propUp = 1
  } else {
    propDown = 1
    propUp = lengthUp / lengthDown
  }
  # cols = pal(101) # keep an extra color in the middle for 0
  numDown = round(nPerSide * propDown)
  numUp = round(nPerSide * propUp)
  cols = c(cols[(nPerSide + 1-numDown):nPerSide], cols[nPerSide + 1], cols[(nPerSide + 2):(nPerSide + 1 + numUp)])
  
  # finally, plot the relative differences
  expectedTypeText = " full debiased"
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    
    if(includeCluster) {
      # gather and plot the debiased results
      nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
      out = load(paste0(nameRoot, '.RData'))
      
      urbanText = ifelse(includeUrban, "", "noUrb")
      clusterText = ifelse(includeCluster, "", "NoClust")
      neither = !includeUrban && !includeUrban
      if(neither)
        debiasedText = "Debiased"
      else
        debiasedText = "debiased"
      typeText = paste0(" ", urbanText, clusterText, debiasedText)
      
      if(i == 1) {
        expected = expit(designRes$predictions$mean)
      }
      vals = expit(designRes$predictions$mean)
      plotMapDat(adm1, plotVar=(vals - expected) / expected, new = TRUE, main=paste0("BYM2 ", typeText, " ", varName, " relative to full debiased"), 
                 cols=cols, zlim=zlim)
    }
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    vals = expit(designRes$predictions$mean)
    plotMapDat(adm1, plotVar=(vals - expected) / expected, new = TRUE, main=paste0("BYM2 ", typeText, " ", varName, " relative to full debiased"), 
               zlim=zlim, cols=cols)
  }
  dev.off()
  
  ## Plot 7: plot the SPDE model relative differences versus the full SPDE model on continuous scale
  # 2 x 2 plot. Columns: cluster, no cluster. Rows: urban, no urban
  argList = list(list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE))
  png(file=paste0("figures/", resultNameRoot, "/fullRelativeSPDE", resultNameRoot, "Continuous.png"), width=800, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
  
  # before the plotting begins, determine the range of relative differences
  zlim=c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    if(i == 1)
      expected = spdeResults$resultsPixel$pred
    
    vals = (spdeResults$resultsPixel$pred - expected) / expected
    zlim = range(c(zlim, vals))
  }
  
  # now shift the color scales so that it center occurs at 0
  cols = relativeCols
  lengthUp = zlim[2] - 0
  lengthDown = 0 - zlim[1]
  if(lengthUp >= lengthDown) {
    # shorten the blue part of the color scale
    propDown = lengthDown / lengthUp
    propUp = 1
  } else {
    propDown = 1
    propUp = lengthUp / lengthDown
  }
  # cols = pal(101) # keep an extra color in the middle for 0
  numDown = round(nPerSide * propDown)
  numUp = round(nPerSide * propUp)
  cols = c(cols[(nPerSide + 1-numDown):nPerSide], cols[nPerSide + 1], cols[(nPerSide + 2):(nPerSide + 1 + numUp)])
  
  # finally, plot the relative differences
  expectedTypeText = " full"
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    vals = spdeResults$resultsPixel$pred
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " ", varName, " relative to full SPDE"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
    quilt.plot(cbind(popGrid$lon, popGrid$lat), (vals-expected)/expected, 
               nx=150, ny=150, add.legend=TRUE, add=TRUE, col=cols, zlim=zlim)
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    # image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
    #            col=cols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 0)
  }
  dev.off()
  
  ## Plot 8: plot the SPDE model relative differences versus the full SPDE model on discrete scale
  # 2 x 2 plot. Columns: cluster, no cluster. Rows: urban, no urban
  argList = list(list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE))
  png(file=paste0("figures/", resultNameRoot, "/fullRelativeSPDE", resultNameRoot, "Discrete.png"), width=800, height=1000)
  par(mfrow=c(2,2))
  
  # before the plotting begins, determine the range of relative differences
  zlim=c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    if(i == 1)
      expected = spdeResults$resultsCounty$pred
    
    vals = (spdeResults$resultsCounty$pred - expected) / expected
    zlim = range(c(zlim, vals))
  }
  
  # now shift the color scales so that it center occurs at 0
  cols = relativeCols
  lengthUp = zlim[2] - 0
  lengthDown = 0 - zlim[1]
  if(lengthUp >= lengthDown) {
    # shorten the blue part of the color scale
    propDown = lengthDown / lengthUp
    propUp = 1
  } else {
    propDown = 1
    propUp = lengthUp / lengthDown
  }
  # cols = pal(101) # keep an extra color in the middle for 0
  numDown = round(nPerSide * propDown)
  numUp = round(nPerSide * propUp)
  cols = c(cols[(nPerSide + 1-numDown):nPerSide], cols[nPerSide + 1], cols[(nPerSide + 2):(nPerSide + 1 + numUp)])
  
  # finally, plot the relative differences
  expectedTypeText = " full"
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    vals = spdeResults$resultsCounty$pred
    plotMapDat(adm1, plotVar=(vals-expected)/expected, new = TRUE, 
               main=paste0("SPDE ", typeText, " ", varName, " relative to full SPDE"), 
               cols=cols, zlim=zlim)
  }
  dev.off()
  
  ## Plot 9: make pair plots of the models
  pdf(file=paste0("figures/", resultNameRoot, "/pairPlot", resultNameRoot, ".pdf"), width=6, height=6)
  # first load full BYM2 and SPDE models
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
  out = load(paste0(nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeUrban
  debiasedText = "Debiased"
  typeTextBYM = paste0(" ", urbanText, clusterText, debiasedText)
  
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeUrban
  notBothText = ifelse(both, "", " ")
  typeTextSPDE = paste0(notBothText, urbanText, clusterText)
  
  my_line <- function(x,y,...){
    points(x,y,..., col="blue")
    abline(a = 0,b = 1,...)
  }
  
  pairs(cbind(naiveResults$est, directEstResults$est, 
              mercerResults$est.mercer, 
              expit(designRes$predictions$mean), spdeResults$resultsCounty$pred), 
        c("Naive", "Direct", "Mercer et al.", "Full BYM2", "Full SPDE"), 
        pch=19, cex=.3, lower.panel=my_line, upper.panel = my_line, 
        main=paste0("County ", varName, " estimate comparisons"))
  dev.off()
  
  ## Plot 10: make pair plots of the BYM2 models
  pdf(file=paste0("figures/", resultNameRoot, "/pairPlotBYM2", resultNameRoot, ".pdf"), width=7, height=7)
  argList = list(list(includeUrbanRural = FALSE, includeCluster = FALSE), 
                 list(includeUrbanRural = FALSE, includeCluster = TRUE), 
                 list(includeUrbanRural = TRUE, includeCluster = FALSE), 
                 list(includeUrbanRural = TRUE, includeCluster = TRUE))
  
  # collect the estimates and model labels
  valMat = c()
  labels = c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    valMat = cbind(valMat, expit(designRes$predictions$mean))
    labels = c(labels, paste0("BYM2", typeText))
    
    if(includeCluster) {
      # also gather and plot the debiased results
      nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
      out = load(paste0(nameRoot, '.RData'))
      
      urbanText = ifelse(includeUrban, "", "noUrb")
      clusterText = ifelse(includeCluster, "", "NoClust")
      both = includeUrban && includeUrban
      if(both)
        debiasedText = "Debiased"
      else
        debiasedText = "debiased"
      typeText = paste0(" ", urbanText, clusterText, debiasedText)
      
      valMat = cbind(valMat, expit(designRes$predictions$mean))
      labels = c(labels, paste0("BYM2", typeText))
    }
  }
  
  # now construct the pair plot
  my_line <- function(x,y,...){
    points(x,y,..., col="blue")
    abline(a = 0,b = 1,...)
  }
  
  pairs(valMat, labels, 
        pch=19, cex=.3, lower.panel=my_line, upper.panel = my_line, 
        main=paste0("BYM2 ", varName, " estimate comparisons"))
  dev.off()
  
  ## Plot 11: make pair plots of the SPDE models region estimates
  pdf(file=paste0("figures/", resultNameRoot, "/pairPlotSPDE", resultNameRoot, "Region.pdf"), width=6, height=6)
  argList = list(list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  # collect the estimates and model labels
  valMat = c()
  labels = c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    vals = spdeResults$resultsRegion$pred
    
    valMat = cbind(valMat, vals)
    labels = c(labels, paste0("SPDE", typeText))
  }
  
  # now construct the pair plot
  my_line <- function(x,y,...){
    points(x,y,..., col="blue")
    abline(a = 0,b = 1,...)
  }
  
  pairs(valMat, labels, 
        pch=19, cex=.8, lower.panel=my_line, upper.panel = my_line, 
        main=paste0("SPDE region ", varName, " estimate comparisons"))
  dev.off()
  
  ## Plot 12: make pair plots of the SPDE models county estimates
  pdf(file=paste0("figures/", resultNameRoot, "/pairPlotSPDE", resultNameRoot, "County.pdf"), width=6, height=6)
  argList = list(list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  # collect the estimates and model labels
  valMat = c()
  labels = c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    vals = spdeResults$resultsCounty$pred
    
    valMat = cbind(valMat, vals)
    labels = c(labels, paste0("SPDE", typeText))
  }
  
  # now construct the pair plot
  my_line <- function(x,y,...){
    points(x,y,..., col="blue")
    abline(a = 0,b = 1,...)
  }
  
  pairs(valMat, labels, 
        pch=19, cex=.3, lower.panel=my_line, upper.panel = my_line, 
        main=paste0("SPDE county ", varName, " estimate comparisons"))
  dev.off()
  
  ## Plot 13: make pair plots of the SPDE models pixel estimates
  pdf(file=paste0("figures/", resultNameRoot, "/pairPlotSPDE", resultNameRoot, "Pixel.pdf"), width=6, height=6)
  argList = list(list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  # collect the estimates and model labels
  valMat = c()
  labels = c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    vals = spdeResults$resultsPixel$pred
    
    valMat = cbind(valMat, vals)
    labels = c(labels, paste0("SPDE", typeText))
  }
  
  # now construct the pair plot
  urban = popGrid$urban
  my_line <- function(x,y,...){
    points(x[!urban],y[!urban],..., col="green")
    points(x[urban],y[urban],..., col="blue")
    abline(a = 0,b = 1,...)
  }
  
  pairs(valMat, labels, 
        pch=".", lower.panel=my_line, upper.panel = my_line, 
        main=paste0("SPDE pixel ", varName, " estimate comparisons"))
  dev.off()
  
  ##### now print the parameter estimates:
  print("printing parameter estimates...")
  
  # BYM2
  argList = list(list(includeUrbanRural = FALSE, includeCluster = FALSE), 
                 list(includeUrbanRural = FALSE, includeCluster = TRUE), 
                 list(includeUrbanRural = TRUE, includeCluster = FALSE), 
                 list(includeUrbanRural = TRUE, includeCluster = TRUE))
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    # add on interval width
    designRes$parameters$width = designRes$parameters$Q90 - designRes$parameters$Q10
    
    print(nameRoot)
    print(xtable(designRes$parameters, digits=3))
  }
  
  # SPDE
  
  argList = list(list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    notBothText = ifelse(both, "", " ")
    typeText = paste0(notBothText, urbanText, clusterText)
    
    print(nameRoot)
    parameters = spdeResults
    parameters$resultsPixel = NULL
    parameters$resultsCounty = NULL
    parameters$resultsRegion = NULL
    parameters = do.call("rbind", parameters)
    print(xtable(parameters, digit=3))
  }
  
  invisible(NULL)
}

makeRedBlueSequentialColors = function(n) {
  # library("colorspace")
  # pal <-choose_palette()
  sequential_hcl(n, h1=10, h2=-115, c1=100, c2=100, l1=44, l2=59, p1=0, p2=2.3)
}

makeRedBlueDivergingColors = function(n) {
  # library("colorspace")
  # pal <-choose_palette()
  diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9)
}

makeBlueSequentialColors = function(n) {
  # library("colorspace")
  # pal <-choose_palette()
  sequential_hcl(n, h1=260, c1=80, l1=30, l2=90, p1=1.5, rev=TRUE)
}

makeBlueYellowSequentialColors = function(n) {
  # library("colorspace")
  # pal <-choose_palette()
  sequential_hcl(n, h1=300, h2=75, c1=40, c2=95, l1=15, l2=90, p1=1.0, p2=1.1)
}

makeRedGreenDivergingColors = function(n) {
  # library("colorspace")
  # pal <-choose_palette()
  diverging_hcl(n, h1=265, h2=101, c1=100, l1=50, l2=92, p1=0.6, p2=1.5)
}