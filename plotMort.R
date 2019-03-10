# script for plotting predictions for neonatal mortality in Kenya

# before we load any results, make a function for what to do with them
plotResults = function(preds, sds, upper90, lower10) {
  plotMapDat(plotVar=preds)
}

# plot the actual data
pdf(file="figures/EAsUrban.pdf", width=5, height=5)
urban = mort$urban
plot(mort$lon[!urban], mort$lat[!urban], pch=19, col="green", main=TeX("Urban vs. rural clusters"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.2)
points(mort$lon[urban], mort$lat[urban], pch=19, col="blue", cex=.2)
# world(add=TRUE)
plotMapDat(adm1)
dev.off()

pdf(file="figures/empiricalMortality.pdf", width=5, height=5)
quilt.plot(mort$lon, mort$lat, mort$y / mort$n, nx=100, ny=100, ylim=kenyaLatRange, xlim=kenyaLonRange, 
           xlab="Longitude", ylab="Latitude", main=TeX("Empirical neonatal mortality rates"))
# world(add=TRUE)
plotMapDat(adm1)
dev.off()

##### Naive and direct estimates
# save(directEstMort, naiveMort, file="resultsDirectNaiveMort.RData")
out = load("resultsDirectNaiveMort.RData")
plotMapDat(adm1, plotVar=naiveMort$u1m, new = TRUE, main="Naive model estimates")
plotMapDat(adm1, plotVar=sqrt(naiveMort$var.est), new = TRUE, main="Naive model logit predictive SDs")
plotMapDat(adm1, plotVar=expit(naiveMort$upper), new = TRUE, main="Naive model 10th percentile")
plotMapDat(adm1, plotVar=expit(naiveMort$lower), new = TRUE, main="Naive model 90th percentile")
plotMapDat(adm1, plotVar=directEstMort$u1m, new = TRUE, main="Direct model estimates")
plotMapDat(adm1, plotVar=sqrt(directEstMort$var.est), new = TRUE, main="Direct model logit predictive SDs")
plotMapDat(adm1, plotVar=expit(directEstMort$upper), new = TRUE, main="Direct model 10th percentile")
plotMapDat(adm1, plotVar=expit(directEstMort$lower), new = TRUE, main="Direct model 90th percentile")

##### Mercer et al. estimates
# save(mercerMort, file=paste0("resultsMercerMort.RData"))
out = load("resultsMercerMort.RData")
plotMapDat(adm1, plotVar=mercerMort$u1m.mercer, new = TRUE, main="Mercer et al. model estimates")
plotMapDat(adm1, plotVar=sqrt(mercerMort$var.est.mercer), new = TRUE, main="Mercer et al. model logit predictive SDs")
plotMapDat(adm1, plotVar=expit(mercerMort$lower.mercer), new = TRUE, main="Mercer et al. model 10th percentile")
plotMapDat(adm1, plotVar=expit(mercerMort$upper.mercer), new = TRUE, main="Mercer et al. model 90th percentile")

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
  
  plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 model estimates", typeText))
  plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 model logit predictive SDs", typeText))
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 model 10th percentile", typeText))
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 model 90th percentile", typeText))
  
  if(includeCluster) {
    # also plot the debiased results
  }
}

##### SPDE estimates
argList = list(list(clustDat = mort, includeClustEffect = FALSE, urbanEffect = FALSE), 
               list(clustDat = mort, includeClustEffect = FALSE, urbanEffect = TRUE), 
               list(clustDat = mort, includeClustEffect = TRUE, urbanEffect = FALSE), 
               list(clustDat = mort, includeClustEffect = TRUE, urbanEffect = TRUE))