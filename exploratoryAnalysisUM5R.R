# Kenya UM5R

setwd("~/Google Drive/UW/Wakefield/WakefieldShared/UM5R/")


# read in the STATA data:
# DHS recode manual:
# https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf

library(haven)
library(fields)

##### first plot our actual dataset

pdf(file="figures/mortDat.pdf", width=5, height=5)
par(mfrow=c(1,1))
quilt.plot(mort$lon, mort$lat, mort$y/mort$n, main=TeX("Survey Empirical UM5R (2003-2007)"), 
           xlab="Longitude", ylab="Latitude", xlim=kenyaLonRange, ylim=kenyaLatRange, nx=120, ny=120)
world(add=TRUE)
plotMapDat(adm1)
dev.off()

pdf(file="figures/EAsUrban.pdf", width=5, height=5)
urban = kenyaEAs$urban
plot(kenyaEAs$lon[!urban], kenyaEAs$lat[!urban], pch=".", col="green", main=TeX("Urban vs. rural enumeration areas"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude")
points(kenyaEAs$lon[urban], kenyaEAs$lat[urban], pch=".", col="blue")
world(add=TRUE)
plotMapDat(adm1)
dev.off()

# lines represent births
data <- data.frame(read_dta("Kenya2014BirthRecode/KEBR70FL.DTA"))

length(unique(data$v001)) # 1593 unqiue clusters
urbRur = by(as.numeric(data$v025), factor(data$v001), mean)
sum(urbRur==1)
# [1] 617
sum(urbRur==2)
# [1] 976

# load Kenya data
out = load("kenyaData.RData")
rm(y)
attach(mort)
out = load("gpsDat.RData")
# out

library(glmm)
mod = glmm(y ~ region + urban, list(~region, ~clusterID), c("region", "clustID"), family=binomial.glmm)

# test urbanicity Data set
library(raster)
library(rgdal)
urban = raster("GHS_SMOD_POP2015_GLOBE_R2016A_54009_1k_v1_0/GHS_SMOD_POP2015_GLOBE_R2016A_54009_1k_v1_0.tif")

# first project to be on same scale as pop dataset
urban = projectRaster(urban, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", method="ngb")
lonRange=c(-180, 180)
latRange=c(-90,90)
kenyaLonRange = c(33.5, 42)
kenyaLatRange = c(-5,5.5)
kenyaLonLength = kenyaLonRange[2] - kenyaLonRange[1]
kenyaLatLength = kenyaLatRange[2] - kenyaLatRange[1]
kenyaExtent = extent(c(xmin=kenyaLonRange[1], xmax=kenyaLonRange[2],
                       ymin=kenyaLatRange[1], ymax=kenyaLatRange[2]))
numEAs = 96251
#  Number of rows and columns at the original resolution
# >  range(rI)
# [1]  805 1008
# >  range(cI)
# [1] 4082 4560
totalRows = 15236 # latitude
totalCols = 35497 # longitude
origLonRes = 360/totalCols # ~98.60/degree
origLatRes = 180/totalRows # ~84.64/degree
resPerDegLon = 1/origLonRes
resPerDegLat = 1/origLatRes
extentCols = round(kenyaLonLength*resPerDegLon)
extentRows = round(kenyaLatLength*resPerDegLat)
increaseFac = 2
lonsInterp = seq(kenyaLonRange[1], kenyaLonRange[2], l=round(extentCols*increaseFac))
latsInterp = seq(kenyaLatRange[1], kenyaLatRange[2], l=round(extentRows*increaseFac))
locsInterp = make.surface.grid(list(x=lonsInterp, y=latsInterp))
# test = extract(urban, kenyaExtent)
kenyaUrbanVals = extract(urban, SpatialPoints(locsInterp),method="simple") # all of these methods have
# kenyaUrbanVals2 = extract(urban, data.frame(lon=locsInterp[,1], lat=locsInterp[,2]), method="simple")
lonRes = lonsInterp[2] - lonsInterp[1]
latRes = latsInterp[2] - latsInterp[1]

sum(!is.na(kenyaUrbanVals))
sum(kenyaUrbanVals != round(kenyaUrbanVals)) # should be 0
kenyaUrban = data.frame(list(lon=locsInterp[,1], lat=locsInterp[,2], urban=kenyaUrbanVals))

# read in administrative map areas
# https://stackoverflow.com/questions/17723822/administrative-regions-map-of-a-country-with-ggmap-and-ggplot2
library(ggplot2)
library(rgdal)
library(sp)
# out = load("mapData/KEN_adm1.rds")
# pakistan.adm2.spdf <- get("gadm")
adm1 = readRDS("mapData/KEN_adm1.rds")
plot(adm1)
adm0 = readRDS("mapData/KEN_adm0.rds")
plot(adm0)
names(adm0)

# subset population density to be within Kenya
polys = adm0@polygons
kenyaPoly = polys[[1]]@Polygons[[77]]@coords
plot(kenyaPoly, type="l")
inKenya = in.poly(cbind(kenyaUrban$lon, kenyaUrban$lat), kenyaPoly)
kenyaUrban = kenyaUrban[inKenya,]
dim(kenyaUrban) # make sure we have much more than 96,000 points
# [1] 396516     3 for increaseFac=1
# [1] 1588080       3 for increaseFac=2

# renormalize Kenya population
kenyaUrban$urbanOrig = kenyaUrban$urban

save(kenyaUrban, file="kenyaUrban.RData")

varRange = range(kenyaUrban$urban)
cols = tim.colors(3)
plotVar = kenyaUrban$urban
varRange=range(plotVar)
#############################################
kenyaEAs = simEAs(kenyaPop)
png("figures/EAsAndUrbanicity.png", width=1000, height=600)
par(mfrow=c(1,2))
# set.panel(1,2)
par( oma=c( 0,0,0,5))
plot(kenyaEAs$lon, kenyaEAs$lat, pch=".", col="blue", main=TeX("Enumeration Areas"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
quilt.plot(kenyaUrban[,1:2], kenyaUrban$urban, 
           nx=400, ny=400, main=TeX("Kenya urbanicity"), ylim=kenyaLatRange, 
           xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
dev.off()

###################################
numClusters = 1593
clustersPerCounty = aggregate(rep(1, nrow(mort)), list(mort$admin1), sum)
clusters = kenyaEAs[sample(1:numEAs, numClusters, replace= FALSE),]
clusters2 = kenyaEAs[sampleByStratum(kenyaEAs$admin1, clustersPerCounty$x),]
png("figures/EAsPopClusters.png", width=1000, height=1000)
par(mfrow=c(2,2))
# set.panel(1,2) I
par( oma=c( 0,0,0,5))
plot(clusters2$lon, clusters2$lat, typ="n", main=TeX("Simulated Clusters (by County)"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
plotMapDat(adm1, lwd=.5)
points(clusters2$lon, clusters2$lat, pch=".", col="blue")
plot(mort$lon, mort$lat, type="n", main=TeX("True Clusters"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
plotMapDat(adm1, lwd=.5)
points(mort$lon, mort$lat, pch=".", col="blue")
plot(clusters$lon, clusters$lat, type="n", main=TeX("Simulated Clusters"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
plotMapDat(adm1, lwd=.5)
points(clusters$lon, clusters$lat, pch=".", col="blue")
plot(kenyaPop[,1:2], type="n", main=TeX("Kenya Population Density (people/mi$^2$)"), ylim=kenyaLatRange, 
     xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
quilt.plot(kenyaPop[,1:2], kenyaPop$pop/cellArea, FUN = function(x) {x = x[x != 0]; log10(mean(x))}, 
           nx=400, ny=400, add.legend=FALSE, add=TRUE)
plotMapDat(adm1, lwd=.5)
world(add=TRUE)
par( oma=c(0,0,0,2))
image.plot(zlim=varRange, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
           col=cols, add = TRUE, axis.args=list(at=log10(ticks), labels=ticks))
dev.off()


### determine the threshold for urban/rural
## first get population density at cluster locations
# clusterPop = extract(pop, SpatialPoints(cbind(mort$lon, mort$lat)),method="bilinear")
clusterUrban = extract(urban, SpatialPoints(cbind(mort$lon, mort$lat)),method="simple")
range(clusterUrban[ mort$urban])
range(clusterUrban[ !mort$urban])
pdf(file="figures/urbanicityByUrbanHist.pdf", width=8, height=5)
par(mfrow=c(1,2))
breaks = seq(-.5, 3.5, by=1)
hist(clusterUrban[ mort$urban], main="Urban Classifications", breaks=breaks, freq=F, xlab="Urbanicity")
hist(clusterUrban[ !mort$urban], main="Rural Classifications", breaks=breaks, freq=F, xlab="Urbanicity")
dev.off()

freqUrban = table(clusterUrban[ mort$urban])
freqRural = table(clusterUrban[ !mort$urban])
total = freqUrban + freqRural
probUrban = freqUrban/total
pdf("figures/empiricalProbUrbanVsUrbanicity.pdf", width=5, height=5)
barplot(probUrban, names.arg=c("0", "1", "2", "3"), main="Empirical probability of urban", 
        xlab="Urbanicity", ylim=c(0,1))
abline(h=.5, lty=2, col="blue")
dev.off()
probUrban
# 0         1         2         3 
# 0.2376238 0.2605042 0.4300000 0.9850746 

png("figures/UrbanicityVsUrban.png", width=1000, height=600)
par(mfrow=c(1,2))
par( oma=c( 0,0,0,5))
plot( mort$lon[ mort$urban], mort$lat[ mort$urban], pch=".",col="blue", ylim=kenyaLatRange, 
      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", main=TeX("DHS Cluster Urban Classification"))
points( mort$lon[ !mort$urban], mort$lat[ !mort$urban], pch=".",col="green")
world(add=TRUE)
quilt.plot(kenyaUrban[,1:2], kenyaUrban$urban, 
           nx=400, ny=400, add.legend=FALSE, main=TeX("Kenya Urbanicity"), ylim=kenyaLatRange, 
           xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
dev.off()
# par( oma=c(0,0,0,2))
# image.plot(zlim=varRange, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
#            col=cols, add = TRUE, axis.args=list(at=log10(ticks), labels=ticks))

par(mfrow=c(1,3))
urbanicity = factor(clusterUrban)
modUrban = mort$urban
modReg = mort$region
modAd = mort$admin1
mod = glm(modUrban ~ clusterUrban, data=mort, family=binomial("probit"))
mod = glm(modUrban ~ clusterUrban + modReg, data=mort, family=binomial("probit"))
mod = glm(modUrban ~ clusterUrban + modAd, data=mort, family=binomial("probit"))
summary( mod)
urbanLevels = c(0, 1, 2, 3)
probs = predict(mod, list(clusterUrban=urbanLevels), type="response")
plot(log10(clusterPop), urban, pch="+")
lines(clustLog10s, probs, col="blue", lwd=2)
thetas = coef(mod)
abline(v=-thetas[1]/thetas[2], col="green") # threshold: logit,probit,cauchit: 4.260,4.269,4.250
sum((modUrban - fitted(mod))^2) # cauchit works the best
sum((modUrban - round(fitted(mod)))^2) # 402 for all of them, or 372 including region, 279 including county


###### plot thresholded population and see if it matches up reasonably with the data
thresh = 18583.71
threshes = c(8000, 10000, 12500, 15000, 17500, 18583)
for(i in 1:length(threshes)) {
  thresh = threshes[i]
  png(paste0("figures/PopThresh", thresh, ".png"), width=1000, height=600)
  par(mfrow=c(1,2))
  par( oma=c( 0,0,0,5))
  plot( mort$lon[ mort$urban], mort$lat[ mort$urban], pch=".",col="blue", ylim=kenyaLatRange, 
        xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", main=TeX("DHS Cluster Urban Classification"))
  points( mort$lon[ !mort$urban], mort$lat[ !mort$urban], pch=".",col="green")
  world(add=TRUE)
  quilt.plot(kenyaPop[,1:2], kenyaPop$popOrig > thresh, FUN=max, col=c("green", "blue"),
             nx=400, ny=400, main=TeX(paste0("Urban Population Threshold (", thresh, ")")), ylim=kenyaLatRange, 
             xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  world(add=TRUE)
  dev.off()
}

##### test sampling clusters by region and urban/rural
par(mfrow=c(1,1))
plot(easpc[,2]/easpc[,4], clustpc[,2]/clustpc[,4], xlim=c(0,1), ylim=c(0,1), col="blue", 
     xlab="Percent urban clusters in county", ylab="Percent urban EAs in county", 
     main="Percent urban EAs vs clusters by county")
abline(0,1)

# threshes = setThresholds2()
# eaThreshes = sapply(1:nrow(kenyaEAs), function(i) {threshes$threshes[threshes$counties == kenyaEAs$admin1[i]]})
# urban = kenyaEAs$popOrig > eaThreshes
pdf(file="figures/EAsClustsUrban.pdf", width=8, height=5)
par(mfrow=c(1,2))
urban = kenyaEAs$urban
plot(kenyaEAs$lon[!urban], kenyaEAs$lat[!urban], pch=".", col="green", main=TeX("Enumeration Areas"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude")
points(kenyaEAs$lon[urban], kenyaEAs$lat[urban], pch=".", col="blue")
world(add=TRUE)
plotMapDat(adm1)

kenyaClusts = simClusters2(kenyaEAs)$clustDat
urban = kenyaClusts$urban
plot(kenyaClusts$lon[!urban], kenyaClusts$lat[!urban], pch=".", col="green", main=TeX("Enumeration Areas"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude")
points(kenyaClusts$lon[urban], kenyaClusts$lat[urban], pch=".", col="blue")
world(add=TRUE)
plotMapDat(adm1)
dev.off()

##### test new thresholds based on population:
threshes = setThresholds2()
clusterPop = extract(pop, SpatialPoints(cbind(mort$lon, mort$lat)),method="bilinear")
urban = setUrbanByThreshes(clusterPop, mort$admin1, threshes)
pdf("figures/urbanicityThreshTest.pdf", width=8, height=5)
par(mfrow=c(1,2))
plot( mort$lon[urban], mort$lat[urban], pch=".",col="blue", ylim=kenyaLatRange, 
      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", main=TeX("Our Cluster Urban Classification"))
points( mort$lon[ !urban], mort$lat[ !urban], pch=".",col="green")
world(add=TRUE)
plotMapDat(adm1)

plot( mort$lon[mort$urban], mort$lat[mort$urban], pch=".",col="blue", ylim=kenyaLatRange, 
      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", main=TeX("DHS Cluster Urban Classification"))
points( mort$lon[ !mort$urban], mort$lat[ !mort$urban], pch=".",col="green")
world(add=TRUE)
plotMapDat(adm1)
dev.off()

##### test simulating data (on logit scale not including nugget should vary between -2.5 and -2)
set.seed(580252)
beta0 = -2
margVar = .15^2
# tausq = .1^2
tausq = .2^2
gamma = -.5
HHoldVar = 0
# HHoldVar = .3^2
HHoldVar = .4^2
urbanOverSample = 2
effRange = 150
out = simDat2(kenyaEAs, clustDat=NULL, nsim=1, urbanOverSample=2,
              beta0=beta0, margVar=margVar, tausq=tausq, gamma=gamma, HHoldVar=HHoldVar, 
              effRange=effRange)
eaDat = out$eaDat
clustDat = out$clustDat

# compute empirical mortality rates
eaPHat = eaDat$died/eaDat$numChildren
clustPHat = clustDat$died/clustDat$numChildren
allDat = c(eaPHat, clustPHat, eaDat$trueProbDeath, eaDat$trueProbDeathNoNug)
allDat = allDat[is.finite(allDat)]
colRange = quantile(allDat, probs=c(0.025, 0.975))

# plot simulated enumeration areas, clusters, and urbanicity
par(mfrow=c(1,2))
plot( eaDat$lon[ !eaDat$urban], eaDat$lat[ !eaDat$urban], pch=".",col="green", ylim=kenyaLatRange, 
      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", main=TeX("EAs and urban classification"))
points( eaDat$lon[eaDat$urban], eaDat$lat[eaDat$urban], pch=".",col="blue")
world(add=TRUE)
plotMapDat(adm1)

plot( clustDat$lon[ !clustDat$urban], clustDat$lat[ !clustDat$urban], pch=".",col="green", ylim=kenyaLatRange, 
      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", main=TeX("Clusters and urban classification"))
points( clustDat$lon[clustDat$urban], clustDat$lat[clustDat$urban], pch=".",col="blue")
world(add=TRUE)
plotMapDat(adm1)

# plot true and empirical probability of death in EAs and clusters along with urbanicity
pdf(file="figures/simDatEmpiricalProb2.pdf", width=10, height=10)
par(mfrow=c(2,2))
par( oma=c( 0,0,0,3))

# true probabilities without nugget/cluster effect
quilt.plot(cbind(eaDat$lon, eaDat$lat), eaDat$trueProbDeathNoNug, zlim=colRange, nx=200, ny=200, FUN = min, 
           main="Median true probabilities", xlab="", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# true probabilities with nugget/cluster effect
quilt.plot(cbind(eaDat$lon, eaDat$lat), eaDat$trueProbDeath, zlim=colRange, nx=200, ny=200, FUN = min, 
           main="True mortality rate (cluster/nugget RE)", xlab="", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# empirical probabilities in enumeration areas
quilt.plot(cbind(eaDat$lon, eaDat$lat), eaPHat, zlim=colRange, nx=50, ny=50, 
           main="EA empirical mortality rate", xlab="Longitude", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# empirical probabilities in clusters
quilt.plot(cbind(clustDat$lon, clustDat$lat), clustPHat, zlim=colRange, nx=50, ny=50, 
           main="Cluster empirical mortality rate", xlab="Longitude", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)
dev.off()

##### how much do we oversample each urban clusters in each county
pdf(file="figures/urbanOversampling.pdf", width=5, height=5)
plot(1,1, main="County urban oversampling", col="blue", type="n", 
     xlab="Proportion or urban EAs", ylab="Proportion of urban clusters", xlim=c(0,1), ylim=c(0,1))
abline(0,1)
points(propUrbPerCounty(eaDat), propUrbPerCounty(clustDat), col="blue")
dev.off()

##### test new cluster sampling scheme (simClusters3) and simDat2
overSampDat = simDat2(kenyaEAs, clustDat=NULL, nsim=4, urbanOverSample=2, 
                      beta0=-2, margVar=.15^2, tausq=.1^2, gamma=-.5)
eaDat = overSampDat$eaDat
clustDat1 = overSampDat$clustDat[[1]]
clustDat2 = overSampDat$clustDat[[2]]
clustDat3 = overSampDat$clustDat[[3]]

# plot 4 simulations
# plot true and empirical probability of death in EAs and clusters along with urbanicity
pdf(file="figures/simDatEmpiricalProb3OverUrban.pdf", width=10, height=10)
par(mfrow=c(2,2))
par( oma=c( 0,0,0,3))

## empirical probabilities in enumeration areas
quilt.plot(cbind(eaDat$lon, eaDat$lat), eaPHat, zlim=colRange, nx=50, ny=50, 
           main="EA empirical mortality rate", xlab="", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

## empirical probabilities in clusters
# sim 1
quilt.plot(cbind(clustDat1$lon, clustDat1$lat), clustPHat, zlim=colRange, nx=50, ny=50, 
           main="Cluster simulation 1", xlab="", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# sim 2
quilt.plot(cbind(clustDat2$lon, clustDat2$lat), clustPHat, zlim=colRange, nx=50, ny=50, 
           main="Cluster simulation 2", xlab="Longitude", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# sim 3
quilt.plot(cbind(clustDat3$lon, clustDat3$lat), clustPHat, zlim=colRange, nx=50, ny=50, 
           main="Cluster simulation 3", xlab="Longitude", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)
dev.off()

##### do the same except without oversampling
tmpClustDat = simClusters3(kenyaEAs, urbanOverSample=1, nsim=4)
clustList = genAndreaFormatFromEAIs(overSampDat$eaDat, tmpClustDat$eaIs, tmpClustDat$sampleWeights)
SRSDat = list(eaDat=overSampDat$eaDat, clustDat=clustList)
eaDat = SRSDat$eaDat
clustDat1 = SRSDat$clustDat[[1]]
clustDat2 = SRSDat$clustDat[[2]]
clustDat3 = SRSDat$clustDat[[3]]

# plot 4 simulations
# plot true and empirical probability of death in EAs and clusters along with urbanicity
pdf(file="figures/simDatEmpiricalProb3SRS.pdf", width=10, height=10)
par(mfrow=c(2,2))
par( oma=c( 0,0,0,3))

## empirical probabilities in enumeration areas
quilt.plot(cbind(eaDat$lon, eaDat$lat), eaPHat, zlim=colRange, nx=50, ny=50, 
           main="EA empirical mortality rate", xlab="", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

## empirical probabilities in clusters
# sim 1
quilt.plot(cbind(clustDat1$lon, clustDat1$lat), clustPHat, zlim=colRange, nx=50, ny=50, 
           main="Cluster simulation 1", xlab="", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# sim 2
quilt.plot(cbind(clustDat2$lon, clustDat2$lat), clustPHat, zlim=colRange, nx=50, ny=50, 
           main="Cluster simulation 2", xlab="Longitude", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# sim 3
quilt.plot(cbind(clustDat3$lon, clustDat3$lat), clustPHat, zlim=colRange, nx=50, ny=50, 
           main="Cluster simulation 3", xlab="Longitude", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)
dev.off()


##### high-res versions:

eaDat = overSampDat$eaDat
clustDat1 = overSampDat$clustDat[[1]]
clustDat2 = overSampDat$clustDat[[2]]
clustDat3 = overSampDat$clustDat[[3]]

# plot 4 simulations
# plot true and empirical probability of death in EAs and clusters along with urbanicity
pdf(file="figures/simDatEmpiricalProb3OverUrbanHR.pdf", width=10, height=10)
par(mfrow=c(2,2))
par( oma=c( 0,0,0,3))

## empirical probabilities in enumeration areas
quilt.plot(cbind(eaDat$lon, eaDat$lat), eaPHat, zlim=colRange, nx=200, ny=200, 
           main="EA empirical mortality rate", xlab="", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

## empirical probabilities in clusters
# sim 1
quilt.plot(cbind(clustDat1$lon, clustDat1$lat), clustPHat, zlim=colRange, nx=70, ny=70, 
           main="Cluster simulation 1", xlab="", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# sim 2
quilt.plot(cbind(clustDat2$lon, clustDat2$lat), clustPHat, zlim=colRange, nx=70, ny=70, 
           main="Cluster simulation 2", xlab="Longitude", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# sim 3
quilt.plot(cbind(clustDat3$lon, clustDat3$lat), clustPHat, zlim=colRange, nx=70, ny=70, 
           main="Cluster simulation 3", xlab="Longitude", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)
dev.off()

##### do the same except without oversampling
eaDat = SRSDat$eaDat
clustDat1 = SRSDat$clustDat[[1]]
clustDat2 = SRSDat$clustDat[[2]]
clustDat3 = SRSDat$clustDat[[3]]

# plot 4 simulations
# plot true and empirical probability of death in EAs and clusters along with urbanicity
pdf(file="figures/simDatEmpiricalProb3SRSHR.pdf", width=10, height=10)
par(mfrow=c(2,2))
par( oma=c( 0,0,0,3))

## empirical probabilities in enumeration areas
quilt.plot(cbind(eaDat$lon, eaDat$lat), eaPHat, zlim=colRange, nx=200, ny=200, 
           main="EA empirical mortality rate", xlab="", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

## empirical probabilities in clusters
# sim 1
quilt.plot(cbind(clustDat1$lon, clustDat1$lat), clustPHat, zlim=colRange, nx=70, ny=70, 
           main="Cluster simulation 1", xlab="", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# sim 2
quilt.plot(cbind(clustDat2$lon, clustDat2$lat), clustPHat, zlim=colRange, nx=70, ny=70, 
           main="Cluster simulation 2", xlab="Longitude", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# sim 3
quilt.plot(cbind(clustDat3$lon, clustDat3$lat), clustPHat, zlim=colRange, nx=70, ny=70, 
           main="Cluster simulation 3", xlab="Longitude", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)
dev.off()

##### show urban oversampling versus SRS variation
overSampDat = simDat2(kenyaEAs, clustDat=NULL, nsim=4, urbanOverSample=2, 
                      beta0=-2, margVar=.15^2, tausq=.1^2, gamma=-.5)
eaDat = overSampDat$eaDat
clustDat = overSampDat$clustDat[[1]]
clustDat2 = overSampDat$clustDat[[2]]
clustDat3 = overSampDat$clustDat[[3]]

tmpClustDat = simClusters3(kenyaEAs, urbanOverSample=1, nsim=4)
clustList = genAndreaFormatFromEAIs(overSampDat$eaDat, tmpClustDat$eaIs, tmpClustDat$sampleWeights)
SRSDat = list(eaDat=overSampDat$eaDat, clustDat=clustList)

pdf(file="figures/clusterSimsAll.pdf", width=8, height=12)
par(mfrow=c(3,2))
clustDat = SRSDat$clustDat[[1]]
urban = clustDat$urban
plot(clustDat$lon[!urban], clustDat$lat[!urban], pch=19, col="green", main=TeX("SRS urban vs. rural 1"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.3)
points(clustDat$lon[urban], clustDat$lat[urban], pch=19, col="blue", cex=.3)
world(add=TRUE)
plotMapDat(adm1)

clustDat = overSampDat$clustDat[[1]]
urban = clustDat$urban
plot(clustDat$lon[!urban], clustDat$lat[!urban], pch=19, col="green", main=TeX("Stratified urban vs. rural 1"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.3)
points(clustDat$lon[urban], clustDat$lat[urban], pch=19, col="blue", cex=.3)
world(add=TRUE)
plotMapDat(adm1)

clustDat = SRSDat$clustDat[[2]]
urban = clustDat$urban
plot(clustDat$lon[!urban], clustDat$lat[!urban], pch=19, col="green", main=TeX("SRS urban vs. rural 2"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.3)
points(clustDat$lon[urban], clustDat$lat[urban], pch=19, col="blue", cex=.3)
world(add=TRUE)
plotMapDat(adm1)

clustDat = overSampDat$clustDat[[2]]
urban = clustDat$urban
plot(clustDat$lon[!urban], clustDat$lat[!urban], pch=19, col="green", main=TeX("Stratified urban vs. rural 2"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.3)
points(clustDat$lon[urban], clustDat$lat[urban], pch=19, col="blue", cex=.3)
world(add=TRUE)
plotMapDat(adm1)

clustDat = SRSDat$clustDat[[3]]
urban = clustDat$urban
plot(clustDat$lon[!urban], clustDat$lat[!urban], pch=19, col="green", main=TeX("SRS urban vs. rural 3"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.3)
points(clustDat$lon[urban], clustDat$lat[urban], pch=19, col="blue", cex=.3)
world(add=TRUE)
plotMapDat(adm1)

clustDat = overSampDat$clustDat[[3]]
urban = clustDat$urban
plot(clustDat$lon[!urban], clustDat$lat[!urban], pch=19, col="green", main=TeX("Stratified urban vs. rural 3"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.3)
points(clustDat$lon[urban], clustDat$lat[urban], pch=19, col="blue", cex=.3)
world(add=TRUE)
plotMapDat(adm1)

dev.off()

overSampDat = simDat2(kenyaEAs, clustDat=NULL, nsim=9, urbanOverSample=2, 
                      beta0=-2, margVar=.15^2, tausq=.1^2, gamma=-.5)
eaDat = overSampDat$eaDat
pdf(file="figures/clusterSimsUrbanOver.pdf", width=12, height=12)
par(mfrow=c(3,3))

for(i in 1:9) {
  clustDat = overSampDat$clustDat[[i]]
  urban = clustDat$urban
  plot(clustDat$lon[!urban], clustDat$lat[!urban], pch=19, col="green", main=TeX(paste0("Stratified urban vs. rural ", i)), xlim=kenyaLonRange, 
       ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.3)
  points(clustDat$lon[urban], clustDat$lat[urban], pch=19, col="blue", cex=.3)
  world(add=TRUE)
  plotMapDat(adm1)
}

dev.off()

tmpClustDat = simClusters3(kenyaEAs, urbanOverSample=1, nsim=9)
clustList = genAndreaFormatFromEAIs(overSampDat$eaDat, tmpClustDat$eaIs, tmpClustDat$sampleWeights)
SRSDat = list(eaDat=overSampDat$eaDat, clustDat=clustList)
pdf(file="figures/clusterSimSRS.pdf", width=12, height=12)
par(mfrow=c(3,3))

for(i in 1:9) {
  clustDat = SRSDat$clustDat[[i]]
  urban = clustDat$urban
  plot(clustDat$lon[!urban], clustDat$lat[!urban], pch=19, col="green", main=TeX(paste0("SRS urban vs. rural ", i)), xlim=kenyaLonRange, 
       ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.3)
  points(clustDat$lon[urban], clustDat$lat[urban], pch=19, col="blue", cex=.3)
  world(add=TRUE)
  plotMapDat(adm1)
}
dev.off()