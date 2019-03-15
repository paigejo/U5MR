# script for plotting predictions for neonatal mortality in Kenya

##### before we make any plots, put all of them on the same scale
##### make multiple scales, two for estimates and quantiles, and two
##### for standard deviations. For each type of scale, make one that 
##### includes the naive and direct estimates, on anther that does not
# naive and direct
out = load("resultsDirectNaiveMort.RData")
meanRange = range(c(expit(naiveMort$upper),expit(naiveMort$lower) ))
meanRange = range(c(meanRange, expit(naiveMort$upper),expit(naiveMort$lower)))
zlim = range(c(expit(directEstMort$upper),expit(directEstMort$lower) ))
meanRange = range(c(meanRange, zlim))
sdRange = range(sqrt(naiveMort$var.est))
zlim = range(sqrt(directEstMort$var.est))
sdRange = range(c(sdRange, zlim))

# make scales specialized for the direct and naive models
meanRangeND = meanRange
sdRangeND = sdRange
sdTicksND = pretty(sdRangeND)
sdTickLabelsND = as.character(sdTicksND)
meanTicksND = c(0.005, 0.01, pretty(c(.01, meanRangeND[2]), n=5))
meanTicksND = meanTicksND[-which(meanTicksND == 0)]
meanTickLabelsND = as.character(meanTicksND)

# mercer
out = load("resultsMercerMort.RData")
zlim = range(c(expit(mercerMort$lower.mercer),expit(mercerMort$upper.mercer)))
meanRange = range(c(meanRange, zlim))
sdRange = range(c(sdRange, sqrt(mercerMort$var.est.mercer)))
meanRange2 = zlim
sdRange2 = range(sqrt(mercerMort$var.est.mercer))

# bym2
argList = list(list(includeUrbanRural = FALSE, includeCluster = FALSE), 
               list(includeUrbanRural = FALSE, includeCluster = TRUE), 
               list(includeUrbanRural = TRUE, includeCluster = FALSE), 
               list(includeUrbanRural = TRUE, includeCluster = TRUE))
for(i in 1:length(argList)) {
  args = argList[[i]]
  includeUrban = args$includeUrbanRural
  includeCluster = args$includeCluster
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0('bym2MortUrbRur',includeUrban, 'Cluster', includeCluster)
  out = load(paste0(nameRoot, '.RData'))
  zlim = range(c(expit(designRes$predictions$Q10),expit(designRes$predictions$Q90)))
  meanRange = range(c(meanRange, zlim))
  sdRange = range(c(sdRange, designRes$predictions$stddev))
  meanRange2 = range(c(meanRange2, zlim))
  sdRange2 = range(c(sdRange2, sqrt(mercerMort$var.est.mercer)))
  
  if(includeCluster) {
    # also gather and plot the debiased results
    nameRoot = paste0('bym2MortUrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
    out = load(paste0(nameRoot, '.RData'))
    
    zlim = range(c(expit(designRes$predictions$Q10),expit(designRes$predictions$Q90)))
    meanRange = range(c(meanRange, zlim))
    sdRange = range(c(sdRange, designRes$predictions$stddev))
    meanRange2 = range(c(meanRange2, zlim))
    sdRange2 = range(c(sdRange2, designRes$predictions$stddev))
  }
}

# spde
argList = list(list(clustDat = mort, includeClustEffect = FALSE, urbanEffect = FALSE), 
               list(clustDat = mort, includeClustEffect = FALSE, urbanEffect = TRUE), 
               list(clustDat = mort, includeClustEffect = TRUE, urbanEffect = FALSE), 
               list(clustDat = mort, includeClustEffect = TRUE, urbanEffect = TRUE))
for(i in 1:length(argList)) {
  args = argList[[i]]
  includeUrban = args$urbanEffect
  includeCluster = args$includeClustEffect
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDEmort_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  zlim = range(c(spdeResults$resultsCounty$lower,spdeResults$resultsCounty$upper))
  meanRange = range(c(meanRange, zlim))
  sdRange = range(c(sdRange, spdeResults$resultsCounty$sds))
  meanRange2 = range(c(meanRange2, zlim))
  sdRange2 = range(c(sdRange2, spdeResults$resultsCounty$sds))
  
  # get a range just for the SPDE continuous
  zlim = range(c(spdeResults$resultsPixel$lower,spdeResults$resultsPixel$upper))
  if(i==1) {
    meanRangeSPDE = zlim
    sdRangeSPDE = range(spdeResults$resultsCounty$sds)
  } else {
    meanRangeSPDE = range(c(meanRangeSPDE, zlim))
    sdRangeSPDE = range(c(sdRangeSPDE, spdeResults$resultsPixel$sds))
  }
}

# set plot tick marks to be reasonable on logit and log scales
sdTicks = pretty(c(.1, sdRange[2]))
sdTickLabels = as.character(sdTicks)
sdTickLabels[c(5, 7)] = ""
meanTicks = pretty(c(.01, meanRange[2]), n=10)
meanTickLabels = as.character(meanTicks)
meanTickLabels[c(5, 7, 9, 11, 13)] = ""
sdTicks2 = pretty(sdRange2)
sdTickLabels2 = as.character(sdTicks2)
meanTicks2 = pretty(meanRange2, n=5)
meanTickLabels2 = as.character(meanTicks2)
meanTicksSPDE = pretty(meanRangeSPDE)
sdTicksSPDE = pretty(sdRangeSPDE)
meanTickLabelsSPDE = as.character(meanTicksSPDE)
sdTickLabelsSPDE = as.character(sdTicksSPDE)

# add in a few extra tick marks
meanTicks = c(.005, meanTicks)
meanTickLabels = c("0.005", meanTickLabels)
sdTicks = c(0.005, 0.01, 0.05, sdTicks)
sdTickLabels = c("0.005", "0.01", "0.05", sdTickLabels)
sdTicks2 = c(0.005, 0.01, 0.05, sdTicks2)
sdTickLabels2 = c("0.005", "0.01", "0.05", sdTickLabels2)

# plot the actual data
png(file="figures/EAsUrban.png", width=500, height=500)
par(oma=c( 0,0,0,3), mar=c(5.1, 4.1, 4.1, 6))
urban = mort$urban
plot(mort$lon[!urban], mort$lat[!urban], pch=19, col="green", main=TeX("Urban vs. rural clusters"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.2)
points(mort$lon[urban], mort$lat[urban], pch=19, col="blue", cex=.2)
# world(add=TRUE)
plotMapDat(adm1)
dev.off()

png(file="figures/empiricalMortality.png", width=500, height=500)
par(oma=c( 0,0,0,3), mar=c(5.1, 4.1, 4.1, 6))
quilt.plot(mort$lon, mort$lat, mort$y / mort$n, nx=150, ny=150, ylim=kenyaLatRange, xlim=kenyaLonRange, 
           xlab="Longitude", ylab="Latitude", main=TeX("Empirical neonatal mortality rates"))
# world(add=TRUE)
plotMapDat(adm1)
dev.off()

png(file="figures/empiricalMortalityLogit.png", width=500, height=500)
par(oma=c( 0,0,0,3), mar=c(5.1, 4.1, 4.1, 6))
cols = tim.colors(30)
ticks = pretty(seq(0, max(mort$y / mort$n), l=5))
ticks = logit(ticks[-1])
varRange = expit(range(ticks))
# par( oma=c( 0,0,0,5)) # save some room for the legend
plot(cbind(mort$lon, mort$lat), type="n", main=TeX("Kenya empirical mortality rates"), ylim=kenyaLatRange, 
     xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
quilt.plot(cbind(mort$lon, mort$lat), logit(mort$y / mort$n), 
           nx=100, ny=100, add.legend=FALSE, add=TRUE, zlim=logit(varRange))
plotMapDat(adm1, lwd=.5)
world(add=TRUE)
# par( oma=c(0,0,0,2))
image.plot(zlim=logit(varRange), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
           col=cols, add = TRUE, axis.args=list(at=ticks, labels=expit(ticks)))
dev.off()

png(file="figures/empiricalMortalityDiscrete.png", width=500, height=500)
par(oma=c( 0,0,0,3), mar=c(5.1, 4.1, 4.1, 6))
totals = aggregate(mort$n, list(mort$admin1), FUN=sum)
counts = aggregate(mort$y, list(mort$admin1), FUN=sum)
plotMapDat(adm1, plotVar=counts$x / totals$x, new = TRUE, main="Empirical neonatal mortality rates")
dev.off()

png(file="figures/denominatorsDiscrete.png", width=500, height=500)
par(oma=c( 0,0,0,3), mar=c(5.1, 4.1, 4.1, 6))
plotMapDat(adm1, plotVar=totals$x, new = TRUE, main="Number of children in sample")
dev.off()

png(file="figures/empiricalMortalityDiscreteLogit.png", width=500, height=500)
par(oma=c( 0,0,0,3), mar=c(5.1, 4.1, 4.1, 6))
totals = aggregate(mort$n, list(mort$admin1), FUN=sum)
counts = aggregate(mort$y, list(mort$admin1), FUN=sum)
plotMapDat(adm1, plotVar=counts$x / totals$x, new = TRUE, main="Empirical neonatal mortality rates", 
           zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
dev.off()

png(file="figures/populationDiscrete.png", width=500, height=500)
countyPops = poppc$popTotal
sortI = sort(poppc$County, index.return=TRUE)$ix
countyPops = countyPops[sortI]
par(oma=c( 0,0,0,3), mar=c(5.1, 4.1, 4.1, 6))
plotMapDat(adm1, plotVar=countyPops, new = TRUE, main="Total population")
dev.off()

out = load("data4directMort.RData")
png(file="figures/extendedMortalityDiscrete.png", width=500, height=500)
testDied = aggregate(resMort$died, list(resMort$admin1), FUN=sum)
testN = aggregate(rep(1, nrow(resMort)), list(resMort$admin1), FUN=sum)
plotMapDat(adm1, plotVar=testDied$x / testN$x, new = TRUE, main="Extended neonatal mortality rates")
dev.off()

cols = tim.colors(64)
thisPop = popGrid$popOrig
totalPop = 43*10^6 # from DHS 2014 survey final report page 2
thisPop = totalPop * thisPop / sum(thisPop) / 5^2 # population density per km^2
png(file=paste0("figures/populationDensity.png"), width=500, height=500)
par(oma=c( 0,0,0,3), mar=c(5.1, 4.1, 4.1, 6))
plot(cbind(popGrid$lon, popGrid$lat), type="n", main=TeX("Population density (people/$km^2$)"), ylim=kenyaLatRange, 
     xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
popRange = log(c(1, max(thisPop)))
popTicks = log(pretty(exp(popRange)))
popTicks = c(log(10), log(100), log(1000), popTicks[-1])
tickLabels = as.character(exp(popTicks))
tickLabels[7] = ""
quilt.plot(cbind(popGrid$lon, popGrid$lat, thisPop), nx=150, ny=150, add.legend=FALSE, add=TRUE, 
           zlim=popRange)
plotMapDat(adm1, lwd=.5)
image.plot(zlim=popRange, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
           col=cols, add = TRUE, axis.args=list(at=popTicks, labels=tickLabels), legend.mar = 0)
dev.off()

png(file=paste0("figures/populationDensityDiscrete.png"), width=500, height=500)
test = aggregate(thisPop*5^2, list(popGrid$admin1), FUN=sum)
plotMapDat(adm1, plotVar=test$x, new = TRUE, main="Integrated population density", zlim=c(0, 4*10^6))
dev.off()

# plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE model logit predictive SDs", typeText), ylim=kenyaLatRange, 
#      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
# quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$sds), 
#            nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(log(sdRange2)))
# plotMapDat(adm1, lwd=.5)
# points(mort$lon, mort$lat, pch=".")
# image.plot(zlim=range(log(sdRange2)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
#            col=cols, add = TRUE, axis.args=list(at=log(sdTicks2), labels=sdTickLabels2), legend.mar = 0)

tab = totals
tab$x = counts$x / totals$x
tab

##### Naive and direct estimates
# save(directEstMort, naiveMort, file="resultsDirectNaiveMort.RData")
out = load("resultsDirectNaiveMort.RData")
plotName = "MortPreds"
png(file=paste0("figures/naive", plotName, ".png"), width=1000, height=1200)
par(mfrow=c(2,2))
zlim = range(c(expit(naiveMort$upper),expit(naiveMort$lower) ))
plotMapDat(adm1, plotVar=naiveMort$u1m, new = TRUE, main="Naive model estimates", zlim=meanRange)
plotMapDat(adm1, plotVar=sqrt(naiveMort$var.est), new = TRUE, main="Naive model logit predictive SDs", zlim=sdRange)
plotMapDat(adm1, plotVar=expit(naiveMort$upper), new = TRUE, main="Naive model 10th percentile", zlim=meanRange)
plotMapDat(adm1, plotVar=expit(naiveMort$lower), new = TRUE, main="Naive model 90th percentile", zlim=meanRange)
dev.off()
png(file=paste0("figures/naive", plotName, "Logit.png"), width=1000, height=1200)
par(mfrow=c(2,2), oma=c( 0,0,0,4))
zlim = range(c(expit(naiveMort$upper),expit(naiveMort$lower) ))
plotMapDat(adm1, plotVar=naiveMort$u1m, new = TRUE, main="Naive model estimates", zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
plotMapDat(adm1, plotVar=sqrt(naiveMort$var.est), new = TRUE, main="Naive model logit predictive SDs", zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels)
plotMapDat(adm1, plotVar=expit(naiveMort$upper), new = TRUE, main="Naive model 10th percentile", zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
plotMapDat(adm1, plotVar=expit(naiveMort$lower), new = TRUE, main="Naive model 90th percentile", zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
dev.off()
png(file=paste0("figures/direct", plotName, ".png"), width=800, height=1000)
par(mfrow=c(2,2))
zlim = range(c(expit(directEstMort$upper),expit(directEstMort$lower) ))
plotMapDat(adm1, plotVar=directEstMort$u1m, new = TRUE, main="Direct model estimates", zlim=meanRange)
plotMapDat(adm1, plotVar=sqrt(directEstMort$var.est), new = TRUE, main="Direct model logit predictive SDs", zlim=sdRange)
plotMapDat(adm1, plotVar=expit(directEstMort$upper), new = TRUE, main="Direct model 10th percentile", zlim=meanRange)
plotMapDat(adm1, plotVar=expit(directEstMort$lower), new = TRUE, main="Direct model 90th percentile", zlim=meanRange)
dev.off()
png(file=paste0("figures/direct", plotName, "Logit.png"), width=800, height=1000)
par(mfrow=c(2,2), oma=c( 0,0,0,4))
plotMapDat(adm1, plotVar=directEstMort$u1m, new = TRUE, main="Direct model estimates", zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
plotMapDat(adm1, plotVar=sqrt(directEstMort$var.est), new = TRUE, main="Direct model logit predictive SDs", zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels)
plotMapDat(adm1, plotVar=expit(directEstMort$upper), new = TRUE, main="Direct model 10th percentile", zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
plotMapDat(adm1, plotVar=expit(directEstMort$lower), new = TRUE, main="Direct model 90th percentile", zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
dev.off()

##### Mercer et al. estimates
# save(mercerMort, file=paste0("resultsMercerMort.RData"))
plotName = "mercerMortPreds"
out = load("resultsMercerMort.RData")
zlim = range(c(expit(mercerMort$lower.mercer),expit(mercerMort$upper.mercer)))
png(file=paste0("figures/", plotName, ".png"), width=800, height=1000)
par(mfrow=c(2,2))
plotMapDat(adm1, plotVar=mercerMort$u1m.mercer, new = TRUE, main="Mercer et al. model estimates", zlim=meanRange)
plotMapDat(adm1, plotVar=sqrt(mercerMort$var.est.mercer), new = TRUE, main="Mercer et al. model logit predictive SDs", zlim=sdRange)
plotMapDat(adm1, plotVar=expit(mercerMort$lower.mercer), new = TRUE, main="Mercer et al. model 10th percentile", zlim=meanRange)
plotMapDat(adm1, plotVar=expit(mercerMort$upper.mercer), new = TRUE, main="Mercer et al. model 90th percentile", zlim=meanRange)
dev.off()

png(file=paste0("figures/", plotName, "Logit.png"), width=800, height=1000)
par(mfrow=c(2,2), oma=c( 0,0,0,4))
plotMapDat(adm1, plotVar=mercerMort$u1m.mercer, new = TRUE, main="Mercer et al. model estimates", zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
plotMapDat(adm1, plotVar=sqrt(mercerMort$var.est.mercer), new = TRUE, main="Mercer et al. model logit predictive SDs", zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels)
plotMapDat(adm1, plotVar=expit(mercerMort$lower.mercer), new = TRUE, main="Mercer et al. model 10th percentile", zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
plotMapDat(adm1, plotVar=expit(mercerMort$upper.mercer), new = TRUE, main="Mercer et al. model 90th percentile", zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
dev.off()

zlim = range(c(expit(mercerMort$lower.mercer),expit(mercerMort$upper.mercer)))
png(file=paste0("figures/", plotName, "2.png"), width=800, height=1000)
par(mfrow=c(2,2))
plotMapDat(adm1, plotVar=mercerMort$u1m.mercer, new = TRUE, main="Mercer et al. model estimates", zlim=meanRange2)
plotMapDat(adm1, plotVar=sqrt(mercerMort$var.est.mercer), new = TRUE, main="Mercer et al. model logit predictive SDs", zlim=sdRange2)
plotMapDat(adm1, plotVar=expit(mercerMort$lower.mercer), new = TRUE, main="Mercer et al. model 10th percentile", zlim=meanRange2)
plotMapDat(adm1, plotVar=expit(mercerMort$upper.mercer), new = TRUE, main="Mercer et al. model 90th percentile", zlim=meanRange2)
dev.off()

png(file=paste0("figures/", plotName, "Logit2.png"), width=800, height=1000)
par(mfrow=c(2,2), oma=c( 0,0,0,4))
plotMapDat(adm1, plotVar=mercerMort$u1m.mercer, new = TRUE, main="Mercer et al. model estimates", zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
plotMapDat(adm1, plotVar=sqrt(mercerMort$var.est.mercer), new = TRUE, main="Mercer et al. model logit predictive SDs", zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2)
plotMapDat(adm1, plotVar=expit(mercerMort$lower.mercer), new = TRUE, main="Mercer et al. model 10th percentile", zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
plotMapDat(adm1, plotVar=expit(mercerMort$upper.mercer), new = TRUE, main="Mercer et al. model 90th percentile", zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
dev.off()

##### BYM2 estimates
# save(file = paste0('bym2MortUrbRur',includeUrbanRural, 'Cluster', includeCluster, '.RData'), 
#      designRes = designRes)
# if(includeCluster) {
#   save(file = paste0('bym2MortUrbRur',includeUrbanRural, 'Cluster', includeCluster, 'debiased.RData'), 
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
  
  nameRoot = paste0('bym2MortUrbRur',includeUrban, 'Cluster', includeCluster)
  out = load(paste0(nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeUrban
  notBothText = ifelse(both, "", " ")
  typeText = paste0(notBothText, urbanText, clusterText)
  
  zlim = range(c(expit(designRes$predictions$Q10),expit(designRes$predictions$Q90)))
  png(file=paste0("figures/", nameRoot, ".png"), width=800, height=1000)
  par(mfrow=c(2,2))
  plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 model estimates", typeText), zlim=meanRange)
  plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 model logit predictive SDs", typeText), zlim=sdRange)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 model 10th percentile", typeText), zlim=meanRange)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 model 90th percentile", typeText), zlim=meanRange)
  dev.off()
  
  png(file=paste0("figures/", nameRoot, "Logit.png"), width=800, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,0,4))
  plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 model estimates", typeText), zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
  plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 model logit predictive SDs", typeText), zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 model 10th percentile", typeText), zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 model 90th percentile", typeText), zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
  dev.off()
  
  png(file=paste0("figures/", nameRoot, "2.png"), width=800, height=1000)
  par(mfrow=c(2,2))
  plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 model estimates", typeText), zlim=meanRange2)
  plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 model logit predictive SDs", typeText), zlim=sdRange2)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 model 10th percentile", typeText), zlim=meanRange2)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 model 90th percentile", typeText), zlim=meanRange2)
  dev.off()
  
  png(file=paste0("figures/", nameRoot, "Logit2.png"), width=800, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,0,4))
  plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 model estimates", typeText), zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
  plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 model logit predictive SDs", typeText), zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 model 10th percentile", typeText), zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 model 90th percentile", typeText), zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
  dev.off()
  
  if(includeCluster) {
    # also gather and plot the debiased results
    nameRoot = paste0('bym2MortUrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
    out = load(paste0(nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeUrban
    if(both)
      debiasedText = "Debiased"
    else
      debiasedText = "debiased"
    typeText = paste0(" ", urbanText, clusterText, debiasedText)
    
    png(file=paste0("figures/", nameRoot, ".png"), width=800, height=1000)
    par(mfrow=c(2,2))
    plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 model estimates", typeText), zlim=meanRange)
    plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 model logit predictive SDs", typeText), zlim=sdRange)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 model 10th percentile", typeText), zlim=meanRange)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 model 90th percentile", typeText), zlim=meanRange)
    dev.off()
    
    png(file=paste0("figures/", nameRoot, "Logit.png"), width=800, height=1000)
    par(mfrow=c(2,2), oma=c( 0,0,0,4))
    plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 model estimates", typeText), zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
    plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 model logit predictive SDs", typeText), zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 model 10th percentile", typeText), zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 model 90th percentile", typeText), zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
    dev.off()
    
    png(file=paste0("figures/", nameRoot, "2.png"), width=800, height=1000)
    par(mfrow=c(2,2))
    plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 model estimates", typeText), zlim=meanRange2)
    plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 model logit predictive SDs", typeText), zlim=sdRange2)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 model 10th percentile", typeText), zlim=meanRange2)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 model 90th percentile", typeText), zlim=meanRange2)
    dev.off()
    
    png(file=paste0("figures/", nameRoot, "Logit2.png"), width=800, height=1000)
    par(mfrow=c(2,2), oma=c( 0,0,0,4))
    plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 model estimates", typeText), zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
    plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 model logit predictive SDs", typeText), zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 model 10th percentile", typeText), zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 model 90th percentile", typeText), zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
    dev.off()
  }
}

##### SPDE estimates
argList = list(list(clustDat = mort, includeClustEffect = FALSE, urbanEffect = FALSE), 
               list(clustDat = mort, includeClustEffect = FALSE, urbanEffect = TRUE), 
               list(clustDat = mort, includeClustEffect = TRUE, urbanEffect = FALSE), 
               list(clustDat = mort, includeClustEffect = TRUE, urbanEffect = TRUE))
for(i in 1:length(argList)) {
  args = argList[[i]]
  includeUrban = args$urbanEffect
  includeCluster = args$includeClustEffect
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDEmort_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeUrban
  notBothText = ifelse(both, "", " ")
  typeText = paste0(notBothText, urbanText, clusterText)
  
  zlim = range(c(spdeResults$resultsCounty$lower,spdeResults$resultsCounty$upper))
  
  png(file=paste0("figures/preds", nameRoot, ".png"), width=800, height=1000)
  par(mfrow=c(2,2))
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE model estimates", typeText), zlim=meanRange)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE model logit predictive SDs", typeText), zlim=sdRange)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE model 10th percentile", typeText), zlim=meanRange)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE model 90th percentile", typeText), zlim=meanRange)
  dev.off()
  
  png(file=paste0("figures/preds", nameRoot, "Logit.png"), width=800, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,0,4))
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE model estimates", typeText), zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE model logit predictive SDs", typeText), zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE model 10th percentile", typeText), zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE model 90th percentile", typeText), zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
  dev.off()
  
  png(file=paste0("figures/preds", nameRoot, "2.png"), width=800, height=1000)
  par(mfrow=c(2,2))
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE model estimates", typeText), zlim=meanRange2)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE model logit predictive SDs", typeText), zlim=sdRange2)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE model 10th percentile", typeText), zlim=meanRange2)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE model 90th percentile", typeText), zlim=meanRange2)
  dev.off()
  
  png(file=paste0("figures/preds", nameRoot, "Logit2.png"), width=800, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,0,4))
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE model estimates", typeText), zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE model logit predictive SDs", typeText), zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE model 10th percentile", typeText), zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE model 90th percentile", typeText), zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
  dev.off()
  
  ## plot continuous prediction surface
  png(file=paste0("figures/preds", nameRoot, ".png"), width=800, height=1000)
  par(mfrow=c(2,2))
  quilt.plot(popGrid$lon, popGrid$lat, spdeResults$resultsPixel$pred, nx=200, ny=200, zlim=meanRange)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE model estimates", typeText), zlim=meanRange)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE model logit predictive SDs", typeText), zlim=sdRange)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE model 10th percentile", typeText), zlim=meanRange)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE model 90th percentile", typeText), zlim=meanRange)
  dev.off()
  png(file=paste0("figures/preds", nameRoot, "Logit.png"), width=800, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,0,4))
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE model estimates", typeText), zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE model logit predictive SDs", typeText), zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE model 10th percentile", typeText), zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE model 90th percentile", typeText), zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
  dev.off()
  
  cols = tim.colors(64)
  png(file=paste0("figures/preds", nameRoot, "ContinuousLogit.png"), width=800, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE model estimates", typeText), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$pred), 
             nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(logit(meanRange)))
  plotMapDat(adm1, lwd=.5)
  points(mort$lon, mort$lat, pch=".")
  image.plot(zlim=range(logit(meanRange)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=cols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels), legend.mar = 0)
  
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE model logit predictive SDs", typeText), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$sds), 
             nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(log(sdRange)))
  plotMapDat(adm1, lwd=.5)
  points(mort$lon, mort$lat, pch=".")
  image.plot(zlim=range(log(sdRange)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=cols, add = TRUE, axis.args=list(at=log(sdTicks), labels=sdTickLabels), legend.mar = 0)
  
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE model 10th percentile", typeText), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$lower), 
             nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(logit(meanRange)))
  plotMapDat(adm1, lwd=.5)
  points(mort$lon, mort$lat, pch=".")
  image.plot(zlim=range(logit(meanRange)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=cols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels), legend.mar = 0)
  
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE model 90th percentile", typeText), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$upper), 
             nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(logit(meanRange)))
  plotMapDat(adm1, lwd=.5)
  points(mort$lon, mort$lat, pch=".")
  image.plot(zlim=range(logit(meanRange)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=cols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels), legend.mar = 5)
  dev.off()
  
  png(file=paste0("figures/preds", nameRoot, "ContinuousLogit2.png"), width=800, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE model estimates", typeText), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$pred), 
             nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(logit(meanRange2)))
  plotMapDat(adm1, lwd=.5)
  points(mort$lon, mort$lat, pch=".")
  image.plot(zlim=range(logit(meanRange2)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=cols, add = TRUE, axis.args=list(at=logit(meanTicks2), labels=meanTickLabels2), legend.mar = 0)
  
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE model logit predictive SDs", typeText), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$sds), 
             nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(log(sdRange2)))
  plotMapDat(adm1, lwd=.5)
  points(mort$lon, mort$lat, pch=".")
  image.plot(zlim=range(log(sdRange2)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=cols, add = TRUE, axis.args=list(at=log(sdTicks2), labels=sdTickLabels2), legend.mar = 0)
  
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE model 10th percentile", typeText), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$lower), 
             nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(logit(meanRange2)))
  plotMapDat(adm1, lwd=.5)
  points(mort$lon, mort$lat, pch=".")
  image.plot(zlim=range(logit(meanRange2)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=cols, add = TRUE, axis.args=list(at=logit(meanTicks2), labels=meanTickLabels2), legend.mar = 0)
  
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE model 90th percentile", typeText), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$upper), 
             nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(logit(meanRange2)))
  plotMapDat(adm1, lwd=.5)
  points(mort$lon, mort$lat, pch=".")
  image.plot(zlim=range(logit(meanRange2)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=cols, add = TRUE, axis.args=list(at=logit(meanTicks2), labels=meanTickLabels2), legend.mar = 5)
  dev.off()
  
  png(file=paste0("figures/preds", nameRoot, "ContinuousLogitSelf.png"), width=800, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE model estimates", typeText), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$pred), 
             nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(logit(meanRangeSPDE)))
  plotMapDat(adm1, lwd=.5)
  points(mort$lon, mort$lat, pch=".")
  image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=cols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 0)
  
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE model logit predictive SDs", typeText), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$sds), 
             nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(log(sdRangeSPDE)))
  plotMapDat(adm1, lwd=.5)
  points(mort$lon, mort$lat, pch=".")
  image.plot(zlim=range(log(sdRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=cols, add = TRUE, axis.args=list(at=log(sdTicksSPDE), labels=sdTickLabelsSPDE), legend.mar = 0)
  
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE model 10th percentile", typeText), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$lower), 
             nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(logit(meanRangeSPDE)))
  plotMapDat(adm1, lwd=.5)
  points(mort$lon, mort$lat, pch=".")
  image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=cols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 0)
  
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE model 90th percentile", typeText), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$upper), 
             nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(logit(meanRangeSPDE)))
  plotMapDat(adm1, lwd=.5)
  points(mort$lon, mort$lat, pch=".")
  image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=cols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 5)
  dev.off()
}

##### now put the predictions from each model together in the same plot when relevant
# put naive and direct together with standard deviations and credible intervals
## make two by four plot with these two models and lower, mean, and upper bounds, and SDs
# put values from the smooth models including all effects together:
## mercer, BYM2 with urban and cluster effects, SPDE with urban and cluster effects
## make three by four plot with these three models and lower, mean, and upper bounds, and SDs

## Plot 1: direct and naive models
out = load("resultsDirectNaiveMort.RData")

png(file=paste0("figures/fullDirectNaive.png"), width=800, height=1200)
par(mfrow=c(4,2), oma=c( 0,0,0,4))
plotMapDat(adm1, plotVar=naiveMort$u1m, new = TRUE, main="Naive model estimates", zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND)
plotMapDat(adm1, plotVar=directEstMort$u1m, new = TRUE, main="Direct model estimates", zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND)
plotMapDat(adm1, plotVar=sqrt(naiveMort$var.est), new = TRUE, main="Naive model logit predictive SDs", zlim=log(sdRangeND), scaleFun=log, scaleFunInverse=exp, ticks=sdTicksND, tickLabels=sdTickLabelsND)
plotMapDat(adm1, plotVar=sqrt(directEstMort$var.est), new = TRUE, main="Direct model logit predictive SDs", zlim=log(sdRangeND), scaleFun=log, scaleFunInverse=exp, ticks=sdTicksND, tickLabels=sdTickLabelsND)
plotMapDat(adm1, plotVar=expit(naiveMort$upper), new = TRUE, main="Naive model 10th percentile", zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND)
plotMapDat(adm1, plotVar=expit(directEstMort$upper), new = TRUE, main="Direct model 10th percentile", zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND)
plotMapDat(adm1, plotVar=expit(naiveMort$lower), new = TRUE, main="Naive model 90th percentile", zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND)
plotMapDat(adm1, plotVar=expit(directEstMort$lower), new = TRUE, main="Direct model 90th percentile", zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND)
dev.off()

## Plot 2: mercer, BYM2 with urban and cluster effects, SPDE with urban and cluster effects
png(file=paste0("figures/fullSmoothed.png"), width=1000, height=1200)
par(mfrow=c(4,3), oma=c( 0,0,0,4))
includeUrban = TRUE
includeCluster = TRUE
clusterText = ifelse(includeCluster, "", "NoClust")

nameRoot = paste0('bym2MortUrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
out = load(paste0(nameRoot, '.RData'))

urbanText = ifelse(includeUrban, "", "noUrb")
clusterText = ifelse(includeCluster, "", "NoClust")
both = includeUrban && includeUrban
debiasedText = "Debiased"
typeTextBYM = paste0(" ", urbanText, clusterText, debiasedText)

includeUrban = TRUE
includeCluster = TRUE
clusterText = ifelse(includeCluster, "", "NoClust")

nameRoot = paste0("SPDEmort_includeClustEffect", includeCluster, 
                  "_urbanEffect", includeUrban)
out = load(paste0("results", nameRoot, '.RData'))

urbanText = ifelse(includeUrban, "", "noUrb")
clusterText = ifelse(includeCluster, "", "NoClust")
both = includeUrban && includeUrban
notBothText = ifelse(both, "", " ")
typeTextSPDE = paste0(notBothText, urbanText, clusterText)

plotMapDat(adm1, plotVar=mercerMort$u1m.mercer, new = TRUE, main="Mercer et al. model estimates", zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 model estimates", typeTextBYM), zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE model estimates", typeTextSPDE), zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)

plotMapDat(adm1, plotVar=sqrt(mercerMort$var.est.mercer), new = TRUE, main="Mercer et al. model logit predictive SDs", zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2)
plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 model logit predictive SDs", typeTextBYM), zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2)
plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE model logit predictive SDs", typeTextSPDE), zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2)

plotMapDat(adm1, plotVar=expit(mercerMort$lower.mercer), new = TRUE, main="Mercer et al. model 10th percentile", zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 model 10th percentile", typeTextBYM), zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE model 10th percentile", typeTextSPDE), zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)

plotMapDat(adm1, plotVar=expit(mercerMort$upper.mercer), new = TRUE, main="Mercer et al. model 90th percentile", zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 model 90th percentile", typeTextBYM), zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE model 90th percentile", typeTextSPDE), zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2)
dev.off()

## Plot 3: plot all of the SPDE plots together (4 x 4 plot)
argList = list(list(clustDat = mort, includeClustEffect = FALSE, urbanEffect = FALSE), 
               list(clustDat = mort, includeClustEffect = FALSE, urbanEffect = TRUE), 
               list(clustDat = mort, includeClustEffect = TRUE, urbanEffect = FALSE), 
               list(clustDat = mort, includeClustEffect = TRUE, urbanEffect = TRUE))
cols = tim.colors(64)
png(file=paste0("figures/preds", nameRoot, "ContinuousTogetherSelf.png"), width=1400, height=1400)
par(mfrow=c(4,4), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
for(i in 1:length(argList)) {
  args = argList[[i]]
  includeUrban = args$urbanEffect
  includeCluster = args$includeClustEffect
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDEmort_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeUrban
  notBothText = ifelse(both, "", " ")
  typeText = paste0(notBothText, urbanText, clusterText)
  
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE model estimates", typeText), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$pred), 
             nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(logit(meanRangeSPDE)))
  plotMapDat(adm1, lwd=.5)
  points(mort$lon, mort$lat, pch=".")
  image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=cols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 0)
}

for(i in 1:length(argList)) {
  args = argList[[i]]
  includeUrban = args$urbanEffect
  includeCluster = args$includeClustEffect
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDEmort_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeUrban
  notBothText = ifelse(both, "", " ")
  typeText = paste0(notBothText, urbanText, clusterText)
  
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE model logit predictive SDs", typeText), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$sds), 
             nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(log(sdRangeSPDE)))
  plotMapDat(adm1, lwd=.5)
  points(mort$lon, mort$lat, pch=".")
  image.plot(zlim=range(log(sdRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=cols, add = TRUE, axis.args=list(at=log(sdTicksSPDE), labels=sdTickLabelsSPDE), legend.mar = 0)
}

for(i in 1:length(argList)) {
  args = argList[[i]]
  includeUrban = args$urbanEffect
  includeCluster = args$includeClustEffect
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDEmort_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeUrban
  notBothText = ifelse(both, "", " ")
  typeText = paste0(notBothText, urbanText, clusterText)
  
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE model 10th percentile", typeText), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$lower), 
             nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(logit(meanRangeSPDE)))
  plotMapDat(adm1, lwd=.5)
  points(mort$lon, mort$lat, pch=".")
  image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=cols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 0)
}

for(i in 1:length(argList)) {
  args = argList[[i]]
  includeUrban = args$urbanEffect
  includeCluster = args$includeClustEffect
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDEmort_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeUrban
  notBothText = ifelse(both, "", " ")
  typeText = paste0(notBothText, urbanText, clusterText)
  
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE model 90th percentile", typeText), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$upper), 
             nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(logit(meanRangeSPDE)))
  plotMapDat(adm1, lwd=.5)
  points(mort$lon, mort$lat, pch=".")
  image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=cols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 5)
}
dev.off()

##### plot relative differences between the models
## all together relative to the bym2
# estimates
library("colorspace")
# pal <-choose_palette()
cols=diverging_hcl(101, h1=265, h2=101, c1=100, l1=50, l2=92, p1=0.6, p2=1.5)
png(file=paste0("figures/fullRelative.png"), width=800, height=1000)
par(mfrow=c(2, 2), oma=c( 0,0,0,4))

zlim = range(c((plotVar=naiveMort$u1m - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), 
               (directEstMort$u1m - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), 
               (mercerMort$u1m.mercer - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), 
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
numDown = round(50 * propDown)
numUp = round(50 * propUp)
cols = c(cols[(51-numDown):50], cols[51], cols[52:(51 + numUp)])
plotMapDat(adm1, plotVar=(plotVar=naiveMort$u1m - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), new = TRUE, main="Naive model relative to BYM2 estimates", zlim=zlim, col=cols)
plotMapDat(adm1, plotVar=(directEstMort$u1m - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), new = TRUE, main="Direct model relative to BYM2 estimates", zlim=zlim, col=cols)
plotMapDat(adm1, plotVar=(mercerMort$u1m.mercer - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), new = TRUE, main="Mercer et al. model relative to BYM2 estimates", zlim=zlim, col=cols)
plotMapDat(adm1, plotVar=(spdeResults$resultsCounty$pred - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), new = TRUE, main=paste0("SPDE model relative to BYM2 estimates", typeTextSPDE), zlim=zlim, col=cols)
dev.off()

##### now print the parameter estimates:
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
  
  nameRoot = paste0('bym2MortUrbRur',includeUrban, 'Cluster', includeCluster)
  out = load(paste0(nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeUrban
  notBothText = ifelse(both, "", " ")
  typeText = paste0(notBothText, urbanText, clusterText)
  
  print(nameRoot)
  print(round(designRes$parameters, 4))
}

# SPDE

argList = list(list(clustDat = mort, includeClustEffect = FALSE, urbanEffect = FALSE), 
               list(clustDat = mort, includeClustEffect = FALSE, urbanEffect = TRUE), 
               list(clustDat = mort, includeClustEffect = TRUE, urbanEffect = FALSE), 
               list(clustDat = mort, includeClustEffect = TRUE, urbanEffect = TRUE))
for(i in 1:length(argList)) {
  args = argList[[i]]
  includeUrban = args$urbanEffect
  includeCluster = args$includeClustEffect
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDEmort_includeClustEffect", includeCluster, 
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