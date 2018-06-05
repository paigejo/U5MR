# cluster variance = 0.01
load("resultsDirectNaiveTausq0.01.RData")
load("resultsMercerTausq0.01.RData")
load("KenyaSpatialDesignResultNewTausq0.01.RData")
load("resultsSPDETausq0.01.RData")

# cluster variance = 0 (no cluster effect)
load("resultsDirectNaiveTausq0.RData")
load("resultsMercerTausq0.RData")
load("KenyaSpatialDesignResultNewTausq0.RData")
load("resultsSPDETausq0.RData")

source("scores.R")

if(0){
  # load("simDataMulti.RData")
  # load a different 1 of these depending on whether a cluster effect should be included 
  # in the simulation of the data or not (tausq is the cluster effect variance)
  load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOver2.RData")
  # load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOver2.RData")
  
  # compute truth based on superpopulation
  regions = sort(unique(SRSDat[[1]]$admin1))
  truthbycounty <- rep(NA, 47)
  
  for(i in 1:47){
    super = overSampDat[[1]][overSampDat[[1]]$admin1 == regions[i],]
    truthbycounty[i] <- sum(super$died)/sum(super$numChildren)
  }
  truth = data.frame(admin1=regions, truth=truthbycounty)
  save(truth, file="truthbycounty.RData")
}

load("truthbycounty.RData")

# compute scores
scoresDirectSRS = scoresNaiveSRS = scoresMercerSRS = scoresBYMSRS = scoresSPDESRS = data.frame()
scoresDirectoverSamp = scoresNaiveoverSamp = scoresMerceroverSamp = scoresBYMoverSamp = scoresSPDEoverSamp = data.frame()
# scoresDirectSRS = scoresNaiveSRS = scoresMercerSRS = scoresBYMSRS = data.frame()
# scoresDirectoverSamp = scoresNaiveoverSamp = scoresMerceroverSamp = scoresBYMoverSamp = data.frame()

for(i in 1:5){
  print(i)
  
  allresSRS = merge(truth, directEstSRS[[i]], by="admin1")
  colnames(allresSRS) = c("admin1", "truth", paste(colnames(allresSRS)[3:8], "direct", sep=""))
  allresSRS = merge(allresSRS, naiveSRS[[i]], by="admin1")
  allresSRS = merge(allresSRS, mercerSRS[[i]], by="admin1")
  allresSRS = merge(allresSRS, spdeSRS[[i]], by="admin1")
  
  allresoverSamp = merge(truth, directEstoverSamp[[i]], by="admin1")
  colnames(allresoverSamp) = c("admin1", "truth", paste(colnames(allresoverSamp)[3:8], "direct", sep=""))
  allresoverSamp = merge(allresoverSamp, naiveoverSamp[[i]], by="admin1")
  allresoverSamp = merge(allresoverSamp, merceroverSamp[[i]], by="admin1")
  allresoverSamp = merge(allresoverSamp, spdeOverSamp[[i]], by="admin1")
  
  useLogit = TRUE
  
  # SRS setting 
  my.biasSRSdirect = bias(logit(allresSRS$truth), allresSRS$logit.estdirect, logit=useLogit)
  my.mseSRSdirect = mse(logit(allresSRS$truth), allresSRS$logit.estdirect, logit=useLogit)
  my.dssSRSdirect = dss(logit(allresSRS$truth), allresSRS$logit.estdirect, allresSRS$var.estdirect)
  my.crpsSRSdirect = crpsNormal(logit(allresSRS$truth), allresSRS$logit.estdirect, allresSRS$var.estdirect)
  my.coverageSRSdirect = coverage(logit(allresSRS$truth), allresSRS$upperdirect, allresSRS$lowerdirect, logit=useLogit)
  
  my.biasSRSnaive = bias(logit(allresSRS$truth), allresSRS$logit.est, logit=useLogit)
  my.mseSRSnaive = mse(logit(allresSRS$truth), allresSRS$logit.est, logit=useLogit)
  my.dssSRSnaive = dss(logit(allresSRS$truth), allresSRS$logit.est, allresSRS$var.est)
  my.crpsSRSnaive = crpsNormal(logit(allresSRS$truth), allresSRS$logit.est, allresSRS$var.est)
  my.coverageSRSnaive = coverage(logit(allresSRS$truth), allresSRS$upper, allresSRS$lower, logit=useLogit)
  
  my.biasSRSmercer = bias(logit(allresSRS$truth), allresSRS$logit.est.mercer, logit=useLogit)
  my.mseSRSmercer = mse(logit(allresSRS$truth), allresSRS$logit.est.mercer, logit=useLogit)
  my.dssSRSmercer = dss(logit(allresSRS$truth), allresSRS$logit.est.mercer, allresSRS$var.est.mercer)
  my.crpsSRSmercer = crpsNormal(logit(allresSRS$truth), allresSRS$logit.est.mercer, allresSRS$var.est.mercer)
  my.coverageSRSmercer = coverage(logit(allresSRS$truth), allresSRS$lower.mercer, allresSRS$upper.mercer, logit=useLogit)
 
  my.biasSRSbym = bias(logit(allresSRS$truth), designRes$SRSdat$mean[,i], logit=useLogit)
  my.mseSRSbym = mse(logit(allresSRS$truth), designRes$SRSdat$mean[,i], logit=useLogit)
  my.dssSRSbym = dss(logit(allresSRS$truth), designRes$SRSdat$mean[,i], (designRes$SRSdat$stddev[,i])^2)
  my.crpsSRSbym = crpsNormal(logit(allresSRS$truth), designRes$SRSdat$mean[,i], (designRes$SRSdat$stddev[,i])^2)
  my.coverageSRSbym = coverage(logit(allresSRS$truth), designRes$SRSdat$Q10[,i],designRes$SRSdat$Q90[,i], logit=useLogit)
  
  my.biasSRSspde = bias(logit(allresSRS$truth), allresSRS$logit.est.spde, logit=useLogit)
  my.mseSRSspde = mse(logit(allresSRS$truth), allresSRS$logit.est.spde, logit=useLogit)
  my.dssSRSspde = dss(logit(allresSRS$truth), allresSRS$logit.est.spde, allresSRS$var.est.spde)
  my.crpsSRSspde = crpsNormal(logit(allresSRS$truth), allresSRS$logit.est.spde, allresSRS$var.est.spde)
  my.coverageSRSspde = coverage(logit(allresSRS$truth), allresSRS$lower.spde, allresSRS$upper.spde, logit=useLogit)
  
  scoresDirectSRS <- rbind(scoresDirectSRS, 
                           data.frame(dataset=i, region=allresSRS$admin1, 
                                      bias=my.biasSRSdirect, 
                                      mse=my.mseSRSdirect,
                                      dss=my.dssSRSdirect,
                                      coverage=my.coverageSRSdirect, 
                                      var=mean(allresSRS$var.estdirect),
                                      crps=my.crpsSRSdirect))
  
  scoresNaiveSRS <- rbind(scoresNaiveSRS, 
                          data.frame(dataset=i, region=allresSRS$admin1, 
                                     bias=my.biasSRSnaive, 
                                     mse=my.mseSRSnaive,
                                     dss=my.dssSRSnaive,
                                     coverage=my.coverageSRSnaive,
                                     var=mean(allresSRS$var.est),
                                     crps=my.crpsSRSnaive))
  scoresMercerSRS <- rbind(scoresMercerSRS, 
                          data.frame(dataset=i, region=allresSRS$admin1, 
                                     bias=my.biasSRSmercer, 
                                     mse=my.mseSRSmercer,
                                     dss=my.dssSRSmercer,
                                     coverage=my.coverageSRSmercer, 
                                     var=mean(allresSRS$var.est.mercer),
                                     crps=my.crpsSRSmercer))

  scoresBYMSRS <- rbind(scoresBYMSRS, 
                           data.frame(dataset=i, region=allresSRS$admin1, 
                                      bias=my.biasSRSbym, 
                                      mse=my.mseSRSbym,
                                      dss=my.dssSRSbym,
                                      coverage=my.coverageSRSbym, 
                                      var=mean((designRes$SRSdat$stddev[,i])^2),
                                      crps=my.crpsSRSbym))
  scoresSPDESRS <- rbind(scoresSPDESRS, 
                           data.frame(dataset=i, region=allresSRS$admin1, 
                                      bias=my.biasSRSspde, 
                                      mse=my.mseSRSspde,
                                      dss=my.dssSRSspde,
                                      coverage=my.coverageSRSspde, 
                                      var=mean(allresSRS$var.est.spde),
                                      crps=my.crpsSRSspde))
    
  # oversampling setting
  my.biasoverSampdirect = bias(logit(allresoverSamp$truth), allresoverSamp$logit.estdirect, logit=useLogit)
  my.mseoverSampdirect = mse(logit(allresoverSamp$truth), allresoverSamp$logit.estdirect, logit=useLogit)
  my.dssoverSampdirect = dss(logit(allresoverSamp$truth), allresoverSamp$logit.estdirect, allresoverSamp$var.estdirect)
  my.crpsoverSampdirect = crpsNormal(logit(allresoverSamp$truth), allresoverSamp$logit.estdirect, allresoverSamp$var.estdirect)
  my.coverageoverSampdirect = coverage(logit(allresoverSamp$truth), allresoverSamp$upperdirect, allresoverSamp$lowerdirect, logit=useLogit)
  
  my.biasoverSampnaive = bias(logit(allresoverSamp$truth), allresoverSamp$logit.est, logit=useLogit)
  my.mseoverSampnaive = mse(logit(allresoverSamp$truth), allresoverSamp$logit.est, logit=useLogit)
  my.dssoverSampnaive = dss(logit(allresoverSamp$truth), allresoverSamp$logit.est, allresoverSamp$var.est)
  my.crpsoverSampnaive = crpsNormal(logit(allresoverSamp$truth), allresoverSamp$logit.est, allresoverSamp$var.est)
  my.coverageoverSampnaive = coverage(logit(allresoverSamp$truth), allresoverSamp$upper, allresoverSamp$lower, logit=useLogit)

  my.biasoverSampmercer = bias(logit(allresoverSamp$truth), allresoverSamp$logit.est.mercer, logit=useLogit)
  my.mseoverSampmercer = mse(logit(allresoverSamp$truth), allresoverSamp$logit.est.mercer, logit=useLogit)
  my.dssoverSampmercer = dss(logit(allresoverSamp$truth), allresoverSamp$logit.est.mercer, allresoverSamp$var.est.mercer)
  my.crpsoverSampmercer = crpsNormal(logit(allresoverSamp$truth), allresoverSamp$logit.est.mercer, allresoverSamp$var.est.mercer)
  my.coverageoverSampmercer = coverage(logit(allresoverSamp$truth), allresoverSamp$lower.mercer, allresoverSamp$upper.mercer, logit=useLogit)

  my.biasoverSampbym = bias(logit(allresoverSamp$truth), designRes$overSampDat$mean[,i], logit=useLogit)
  my.mseoverSampbym = mse(logit(allresoverSamp$truth), designRes$overSampDat$mean[,i], logit=useLogit)
  my.dssoverSampbym = dss(logit(allresoverSamp$truth), designRes$overSampDat$mean[,i], (designRes$overSampDat$stddev[,i])^2)
  my.crpsoverSampbym = crpsNormal(logit(allresoverSamp$truth), designRes$overSampDat$mean[,i], (designRes$overSampDat$stddev[,i])^2)
  my.coverageoverSampbym = coverage(logit(allresoverSamp$truth), designRes$overSampDat$Q10[,i],designRes$overSampDat$Q90[,i], logit=useLogit)
  
  my.biasoverSampspde = bias(logit(allresoverSamp$truth), allresoverSamp$logit.est.spde, logit=useLogit)
  my.mseoverSampspde = mse(logit(allresoverSamp$truth), allresoverSamp$logit.est.spde, logit=useLogit)
  my.dssoverSampspde = dss(logit(allresoverSamp$truth), allresoverSamp$logit.est.spde, allresoverSamp$var.est.spde)
  my.crpsoverSampspde = crpsNormal(logit(allresoverSamp$truth), allresoverSamp$logit.est.spde, allresoverSamp$var.est.spde)
  my.coverageoverSampspde = coverage(logit(allresoverSamp$truth), allresoverSamp$lower.spde, allresoverSamp$upper.spde, logit=useLogit)
  
  scoresDirectoverSamp <- rbind(scoresDirectoverSamp, 
                                data.frame(dataset=i, region=allresoverSamp$admin1, 
                                           bias=my.biasoverSampdirect, 
                                           mse=my.mseoverSampdirect,
                                           dss=my.dssoverSampdirect,
                                           coverage=my.coverageoverSampdirect, 
                                           var=mean(allresoverSamp$var.estdirect),
                                           crps=my.crpsoverSampdirect))
  scoresNaiveoverSamp<- rbind(scoresNaiveoverSamp, 
                              data.frame(dataset=i, region=allresoverSamp$admin1, 
                                         bias=my.biasoverSampnaive, 
                                         mse=my.mseoverSampnaive,
                                         dss=my.dssoverSampnaive,
                                         coverage=my.coverageoverSampnaive, 
                                         var=mean(allresoverSamp$var.est),
                                         crps=my.crpsoverSampnaive))
  scoresMerceroverSamp<- rbind(scoresMerceroverSamp, 
                              data.frame(dataset=i, region=allresoverSamp$admin1, 
                                         bias=my.biasoverSampmercer, 
                                         mse=my.mseoverSampmercer,
                                         dss=my.dssoverSampmercer,
                                         coverage=my.coverageoverSampmercer, 
                                         var=mean(allresoverSamp$var.est.mercer),
                                         crps=my.crpsoverSampmercer))

  scoresBYMoverSamp <- rbind(scoresBYMoverSamp, 
                        data.frame(dataset=i, region=allresSRS$admin1, 
                                   bias=my.biasoverSampbym, 
                                   mse=my.mseoverSampbym,
                                   dss=my.dssoverSampbym,
                                   coverage=my.coverageoverSampbym, 
                                   var=mean((designRes$overSampDat$stddev[,i])^2),
                                   crps=my.crpsoverSampbym))
  
  scoresSPDEoverSamp <- rbind(scoresSPDEoverSamp, 
                             data.frame(dataset=i, region=allresSRS$admin1, 
                                        bias=my.biasoverSampspde, 
                                        mse=my.mseoverSampspde,
                                        dss=my.dssoverSampspde,
                                        coverage=my.coverageoverSampspde, 
                                        var=mean(allresoverSamp$var.est.spde),
                                        crps=my.crpsoverSampspde))
}

# compare all four 
# pdf("Figures/biasbyregionoverSamp.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(bias~region, data=scoresDirectoverSamp, at=seq(-1, 229, by=5),
#         col="yellow", xlim=c(0,232), names=FALSE, xaxt="n")
# boxplot(bias~region, data=scoresNaiveoverSamp, at=seq(0, 230, by=5), col="orange", xlim=c(0,230), add=TRUE)
# boxplot(bias~region, data=scoresMerceroverSamp, at=seq(1, 231, by=5), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
# boxplot(bias~region, data=scoresBYMoverSamp, at=seq(2, 232, by=5), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
# axis(2, at=seq(0.5, 230.5, by=5), labels=scoresDirectoverSamp$region[1:47])
# legend("top", c("Direct estimates", "Naive", "Mercer", "BYM"),
#        fill = c("yellow", "orange", "green", "lightblue"), ncol=4, cex=2)
# abline(h=0, lwd=2, col=2)
# dev.off()
# 
# pdf("Figures/biasbyregionSRS.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(bias~region, data=scoresDirectoverSamp, at=seq(-1, 229, by=5), 
#         col="yellow", xlim=c(0,232), names=FALSE, xaxt="n")
# boxplot(bias~region, data=scoresNaiveSRS, at=seq(0, 230, by=5), col="orange", xlim=c(0,230), add=TRUE)
# boxplot(bias~region, data=scoresMercerSRS, at=seq(1, 231, by=5), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
# boxplot(bias~region, data=scoresBYMSRS, at=seq(2, 232, by=5), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
# axis(2, at=seq(0.5, 230.5, by=5), labels=scoresDirectSRS$region[1:47])
# legend("top", c("Direct estimates", "Naive", "Mercer", "BYM"),
#        fill = c("yellow", "orange", "green", "lightblue"), ncol=4, cex=2)
# abline(h=0, lwd=2, col=2)
# dev.off()
# 
# pdf("Figures/crpsbyregionoverSamp.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(crps~region, data=scoresDirectoverSamp, at=seq(-1, 229, by=5), 
#         col="yellow", xlim=c(0,232), names=FALSE, xaxt="n")
# boxplot(crps~region, data=scoresNaiveoverSamp, at=seq(0, 230, by=5), col="orange", xlim=c(0,230), add=TRUE)
# boxplot(crps~region, data=scoresMerceroverSamp, at=seq(1, 231, by=5), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
# boxplot(crps~region, data=scoresBYMoverSamp, at=seq(2, 232, by=5), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
# axis(2, at=seq(0.5, 230.5, by=5), labels=scoresDirectoversamp$region[1:47])
# legend("top", c("Direct estimates", "Naive", "Mercer", "BYM"),
#        fill = c("yellow", "orange", "green", "lightblue"), ncol=4, cex=2)
# dev.off()
# 
# pdf("Figures/crpsbyregionSRS.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(crps~region, data=scoresDirectoverSamp, at=seq(-1, 229, by=5), 
#         col="yellow", xlim=c(0,232), names=FALSE, xaxt="n")
# boxplot(crps~region, data=scoresNaiveSRS, at=seq(0, 230, by=5), col="orange", xlim=c(0,230), add=TRUE)
# boxplot(crps~region, data=scoresMercerSRS, at=seq(1, 231, by=5), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
# boxplot(crps~region, data=scoresBYMSRS, at=seq(2, 232, by=5), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
# axis(2, at=seq(0.5, 230.5, by=5), labels=scoresDirectSRS$region[1:47])
# legend("top", c("Direct estimates", "Naive", "Mercer", "BYM"),
#        fill = c("yellow", "orange", "green", "lightblue"), ncol=4, cex=2)
# dev.off()

# compare all five
pdf("figures/biasbyregionoverSamp.pdf", width=20, height=12)
par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
boxplot(bias~region, data=scoresDirectoverSamp, at=seq(-1, 275, by=6), 
        col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
boxplot(bias~region, data=scoresNaiveoverSamp, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
boxplot(bias~region, data=scoresMerceroverSamp, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
boxplot(bias~region, data=scoresBYMoverSamp, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
boxplot(bias~region, data=scoresSPDEoverSamp, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
# axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirectoverSamp$region[1:47])
# axis(1, at=seq(0.5, 276.5, by=6), labels=scoresDirectoverSamp$region[1:47])
legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
       fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
abline(h=0, lwd=2, col=2)
dev.off()

pdf("figures/biasbyregionSRS.pdf", width=20, height=12)
par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
boxplot(bias~region, data=scoresDirectoverSamp, at=seq(-1, 275, by=6), 
        col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
boxplot(bias~region, data=scoresNaiveSRS, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
boxplot(bias~region, data=scoresMercerSRS, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
boxplot(bias~region, data=scoresBYMSRS, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
boxplot(bias~region, data=scoresSPDESRS, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
# axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirectSRS$region[1:47])
# axis(1, at=seq(0.5, 276.5, by=6), labels=scoresDirectSRS$region[1:47])
legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
       fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
abline(h=0, lwd=2, col=2)
dev.off()

pdf("figures/crpsbyregionoverSamp.pdf", width=20, height=12)
par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
boxplot(crps~region, data=scoresDirectoverSamp, at=seq(-1, 275, by=6), 
        col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
boxplot(crps~region, data=scoresNaiveoverSamp, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
boxplot(crps~region, data=scoresMerceroverSamp, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
boxplot(crps~region, data=scoresBYMoverSamp, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
boxplot(crps~region, data=scoresSPDEoverSamp, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
# axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirectoverSamp$region[1:47])
legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
       fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
dev.off()

pdf("figures/crpsbyregionSRS.pdf", width=20, height=12)
par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
boxplot(crps~region, data=scoresDirectoverSamp, at=seq(-1, 275, by=6), 
        col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
boxplot(crps~region, data=scoresNaiveSRS, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
boxplot(crps~region, data=scoresMercerSRS, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
boxplot(crps~region, data=scoresBYMSRS, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
boxplot(crps~region, data=scoresSPDESRS, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
# axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirectSRS$region[1:47])
legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
       fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
dev.off()

# 
# 
# # compare naive, direct and mercer
# pdf("Figures/biasbyregionoverSamp.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(bias~region, data=scoresDirectoverSamp, at=seq(-1, 183, by=4), 
#         col="yellow", xlim=c(0,188), names=FALSE, xaxt="n")
# boxplot(bias~region, data=scoresNaiveoverSamp, at=seq(0, 184, by=4), col="orange", xlim=c(0,138), add=TRUE)
# boxplot(bias~region, data=scoresMerceroverSamp, at=seq(1, 185, by=4), col="green", xlim=c(0,138), add=TRUE, xaxt="n")
# axis(2, at=seq(0, 184, by=4), labels=scoresDirectoversamp$region[1:47])
# legend("top", c("Direct estimates", "Naive", "Mercer"),
#        fill = c("yellow", "orange", "green"), ncol=3, cex=2)
# abline(h=0, lwd=2, col=2)
# dev.off()
# 
# 
# pdf("Figures/biasbyregionSRS.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(bias~region, data=scoresDirectSRS, at=seq(-1, 183, by=4), 
#         col="yellow", xlim=c(0,188), names=FALSE, xaxt="n")
# boxplot(bias~region, data=scoresNaiveSRS, at=seq(0, 184, by=4), col="orange", xlim=c(0,138), add=TRUE)
# boxplot(bias~region, data=scoresMercerSRS, at=seq(1, 185, by=4), col="green", xlim=c(0,138), add=TRUE, xaxt="n")
# axis(2, at=seq(0, 184, by=4), labels=scoresDirectoversamp$region[1:47])
# legend("top", c("Direct estimates", "Naive", "Mercer"),
#        fill = c("yellow", "orange", "green"), ncol=3, cex=2)
# abline(h=0, lwd=2, col=2)
# dev.off()
# 
# 
# pdf("Figures/crpsbyregionoverSamp.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(crps~region, data=scoresDirectoverSamp, at=seq(-1, 183, by=4), 
#         col="yellow", xlim=c(0,188), names=FALSE, xaxt="n")
# boxplot(crps~region, data=scoresNaiveoverSamp, at=seq(0, 184, by=4), col="orange", xlim=c(0,138), add=TRUE)
# boxplot(crps~region, data=scoresMerceroverSamp, at=seq(1, 185, by=4), col="green", xlim=c(0,138), add=TRUE, xaxt="n")
# axis(2, at=seq(0, 184, by=4), labels=scoresDirectoversamp$region[1:47])
# legend("top", c("Direct estimates", "Naive", "Mercer"),
#        fill = c("yellow", "orange", "green"), ncol=3, cex=2)
# dev.off()
# 
# 
# pdf("Figures/crpsbyregionSRS.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(crps~region, data=scoresDirectSRS, at=seq(-1, 183, by=4), 
#         col="yellow", xlim=c(0,188), names=FALSE, xaxt="n")
# boxplot(crps~region, data=scoresNaiveSRS, at=seq(0, 184, by=4), col="orange", xlim=c(0,138), add=TRUE)
# boxplot(crps~region, data=scoresMercerSRS, at=seq(1, 185, by=4), col="green", xlim=c(0,138), add=TRUE, xaxt="n")
# axis(2, at=seq(0, 184, by=4), labels=scoresDirectoversamp$region[1:47])
# legend("top", c("Direct estimates", "Naive", "Mercer"),
#        fill = c("yellow", "orange", "green"), ncol=3, cex=2)
# dev.off()

# comparing only naive to direct
# 
# pdf("Figures/biasbyregionoverSamp.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(bias~region, data=scoresDirectoverSamp, at=seq(-1, 137, by=3), 
#         col="yellow", xlim=c(0,138), names=FALSE, xaxt="n")
# boxplot(bias~region, data=scoresNaiveoverSamp, at=seq(0, 138, by=3), col="orange", xlim=c(0,138), add=TRUE)
# axis(2, at=seq(-0.5, 137.5, by=3), labels=scoresDirectoversamp$region[1:47])
# legend("top", c("Direct estimates", "Naive"),
#        fill = c("yellow", "orange"), ncol=2, cex=2)
# abline(h=0, lwd=2, col=2)
# dev.off()

# pdf("Figures/biasbyregionSRS.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(bias~region, data=scoresDirectSRS, at=seq(-1, 137, by=3), 
#         col="yellow", xlim=c(0,138), names=FALSE, xaxt="n")
# boxplot(bias~region, data=scoresNaiveSRS, at=seq(0, 138, by=3), col="orange", xlim=c(0,138), add=TRUE)
# axis(2, at=seq(-0.5, 137.5, by=3), labels=scoresDirectSRS$region[1:47])
# legend("top", c("Direct estimates", "Naive"),
#        fill = c("yellow", "orange"), ncol=2, cex=2)
# abline(h=0, lwd=2, col=2)
# dev.off()

# 
# pdf("Figures/crpsbyregionoverSamp.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(crps~region, data=scoresDirectoverSamp, at=seq(-1, 137, by=3), 
#         col="yellow", xlim=c(0,138), names=FALSE, xaxt="n")
# boxplot(crps~region, data=scoresNaiveoverSamp, at=seq(0, 138, by=3), col="orange", xlim=c(0,138), add=TRUE)
# axis(2, at=seq(-0.5, 137.5, by=3), labels=scoresDirectoversamp$region[1:47])
# legend("top", c("Direct estimates", "Naive"),
#        fill = c("yellow", "orange"), ncol=2, cex=2)
# dev.off()
# 
# 
# pdf("Figures/crpsbyregionSRS.pdf", width=20, height=12)
# par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
# boxplot(crps~region, data=scoresDirectSRS, at=seq(-1, 137, by=3), 
#         col="yellow", xlim=c(0,138), names=FALSE, xaxt="n")
# boxplot(crps~region, data=scoresNaiveSRS, at=seq(0, 138, by=3), col="orange", xlim=c(0,138), add=TRUE)
# axis(2, at=seq(-0.5, 137.5, by=3), labels=scoresDirectSRS$region[1:47])
# legend("top", c("Direct estimates", "Naive"),
#        fill = c("yellow", "orange"), ncol=2, cex=2)
# dev.off()

naive = apply(scoresNaiveSRS[, c("bias", "mse", "dss", "crps", "var", "coverage")], 2, mean)
direct = apply(scoresDirectSRS[, c("bias", "mse", "dss", "crps", "var","coverage")], 2, mean)
mercer = apply(scoresMercerSRS[, c("bias", "mse", "dss", "crps", "var","coverage")], 2, mean)
bym = apply(scoresBYMSRS[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean)
spde = apply(scoresSPDESRS[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean)
idx = c(1,2,5,4,6)
tab = rbind(c( naive[idx]),
      c( direct[idx]),
      c( mercer[idx]),
      c( bym[idx]), 
      c( spde[idx]))
rownames(tab) = c("Naive", "Direct estimates", "Mercer et al.", "Model-based BYM", "SPDE")
library(xtable)
print(xtable(tab, digits=3), 
      only.contents=TRUE, 
      include.colnames=FALSE,
      hline.after=NULL)
naive = round(apply(scoresNaiveoverSamp[, c("bias","mse", "dss",  "crps", "var","coverage")], 2, mean),3)
direct = round(apply(scoresDirectoverSamp[, c("bias","mse", "dss",  "crps","var","coverage")], 2, mean),3)
mercer = round(apply(scoresMerceroverSamp[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean),3)
bym = round(apply(scoresBYMoverSamp[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean),3)
spde = round(apply(scoresSPDEoverSamp[, c("bias", "mse", "dss", "crps","var", "coverage")], 2, mean),3)
tab = rbind(c( naive[idx]),
            c( direct[idx]),
            c( mercer[idx]),
            c( bym[idx]), 
            c( spde[idx]))
rownames(tab) = c("Naive", "Direct estimates", "Mercer et al.", "Model-based BYM", "SPDE")
library(xtable)
print(xtable(tab, digits=3), only.contents=TRUE, 
      include.rownames=FALSE,
      include.colnames=FALSE,
      hline.after=NULL)




naive = apply(scoresNaiveSRS[, c("bias", "mse",  "coverage")], 2, mean)
direct = apply(scoresDirectSRS[, c("bias", "mse", "coverage")], 2, mean)
mercer = apply(scoresMercerSRS[, c("bias", "mse", "coverage")], 2, mean)
bym = apply(scoresBYMSRS[, c("bias", "mse", "coverage")], 2, mean)
spde = apply(scoresSPDESRS[, c("bias", "mse", "coverage")], 2, mean)
idx = c(1,2)
tab = rbind(c( naive[idx]),
            c( direct[idx]),
            c( mercer[idx]),
            c( bym[idx]), 
            c( spde[idx]))
rownames(tab) = c("Naive", "Direct estimates", "Mercer et al.", "Model-based BYM", "SPDE")
library(xtable)
print(xtable(tab, digits=3), 
      only.contents=TRUE, 
      include.colnames=TRUE,
      hline.after=NULL)

naive = round(apply(scoresNaiveoverSamp[, c("bias","mse", "coverage")], 2, mean),3)
direct = round(apply(scoresDirectoverSamp[, c("bias","mse", "coverage")], 2, mean),3)
mercer = round(apply(scoresMerceroverSamp[, c("bias", "mse", "coverage")], 2, mean),3)
bym = round(apply(scoresBYMoverSamp[, c("bias", "mse", "coverage")], 2, mean),3)
spde = round(apply(scoresSPDEoverSamp[, c("bias", "mse", "coverage")], 2, mean),3)

tab = rbind(c( naive[idx]),
            c( direct[idx]),
            c( mercer[idx]),
            c( bym[idx]), 
            c( spde[idx]))
rownames(tab) = c("Naive", "Direct estimates", "Mercer et al.", "Model-based BYM", "SPDE")
library(xtable)
print(xtable(tab, digits=3), only.contents=TRUE, 
      include.rownames=TRUE,
      include.colnames=TRUE,
      hline.after=NULL)
