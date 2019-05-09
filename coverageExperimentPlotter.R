# this script plots results from the coverage experiment
library(fields)
library(RColorBrewer)

# first load the experiment results
out = load("coverageSimulation.RData")

## plot coverage versus s*
pdf("figures/coverageExperiment/cov50sStar.pdf", width=5, height=5)
plot(ss, apply(cov50, 4, mean), type="l", main="50% Coverage Versus s", xlab="s", ylab="Coverage Probability", 
     ylim=c(0, 1), col="blue")
abline(h=.5, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov60sStar.pdf", width=5, height=5)
plot(ss, apply(cov60, 4, mean), type="l", main="60% Coverage Versus s", xlab="s", ylab="Coverage Probability", 
     ylim=c(0, 1), col="blue")
abline(h=.6, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov70sStar.pdf", width=5, height=5)
plot(ss, apply(cov70, 4, mean), type="l", main="70% Coverage Versus s", xlab="s", ylab="Coverage Probability", 
     ylim=c(0, 1), col="blue")
abline(h=.7, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov80sStar.pdf", width=5, height=5)
plot(ss, apply(cov80, 4, mean), type="l", main="80% Coverage Versus s", xlab="s", ylab="Coverage Probability", 
     ylim=c(0, 1), col="blue")
abline(h=.8, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov90sStar.pdf", width=5, height=5)
plot(ss, apply(cov90, 4, mean), type="l", main="90% Coverage Versus s", xlab="s", ylab="Coverage Probability", 
     ylim=c(0, 1), col="blue")
abline(h=.9, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov95sStar.pdf", width=5, height=5)
plot(ss, apply(cov95, 4, mean), type="l", main="95% Coverage Versus s", xlab="s", ylab="Coverage Probability", 
     ylim=c(0, 1), col="blue")
abline(h=.95, lty=2)
dev.off()


## plot coverage versus u(s*)
# make the bins over values of u(s*)
# breaks = c(-Inf, seq(min(uSamples), max(uSamples), l=51)[-c(1, 50)], Inf)
breaks = qnorm(seq(1e-06, 1-1e-06, l=40))

# uSamples = array(dim=c(ns, n + m, nu))
# cov50 = array(dim=c(ns, nu, ny, m))
# uHats = array(dim=c(ns, m, nu, ny))
uSamplesPerm = aperm(uSamples, c(1, 3, 2))[,,-(1:n)]
out = stats.bin(uSamplesPerm, cov50y, breaks=breaks)
x = out$centers
y50 = out$stats[2,]
y60 = stats.bin(uSamplesPerm, cov60y, breaks=breaks)$stats[2,]
y70 = stats.bin(uSamplesPerm, cov70y, breaks=breaks)$stats[2,]
y80 = stats.bin(uSamplesPerm, cov80y, breaks=breaks)$stats[2,]
y90 = stats.bin(uSamplesPerm, cov90y, breaks=breaks)$stats[2,]
y95 = stats.bin(uSamplesPerm, cov95y, breaks=breaks)$stats[2,]
pdf("figures/coverageExperiment/cov50u.pdf", width=5, height=5)
plot(x, y50, type="l", main="50% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.5, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov60u.pdf", width=5, height=5)
plot(x, y60, type="l", main="60% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.6, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov70u.pdf", width=5, height=5)
plot(x, y70, type="l", main="70% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.7, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov80u.pdf", width=5, height=5)
plot(x, y80, type="l", main="80% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.8, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov90u.pdf", width=5, height=5)
plot(x, y90, type="l", main="90% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.9, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov95u.pdf", width=5, height=5)
plot(x, y95, type="l", main="95% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.95, lty=2)
dev.off()

## plot coverage versus u(s*), filtered within the data domain
# make the bins over values of u(s*)
# breaks = c(-Inf, seq(min(uSamples), max(uSamples), l=51)[-c(1, 50)], Inf)
breaks = qnorm(seq(1e-06, 1-1e-06, l=50))

# uSamples = array(dim=c(ns, n + m, nu))
# cov50 = array(dim=c(ns, nu, ny, m))
# uHats = array(dim=c(ns, m, nu, ny))
uSamplesPerm = aperm(uSamples, c(1, 3, 2))[,,-(1:n)]
out = stats.bin(uSamplesPerm[,,1:33], cov50y[,,1:33], breaks=breaks)
x = out$centers
y50 = out$stats[2,]
y60 = stats.bin(uSamplesPerm[,,1:33], cov60y[,,1:33], breaks=breaks)$stats[2,]
y70 = stats.bin(uSamplesPerm[,,1:33], cov70y[,,1:33], breaks=breaks)$stats[2,]
y80 = stats.bin(uSamplesPerm[,,1:33], cov80y[,,1:33], breaks=breaks)$stats[2,]
y90 = stats.bin(uSamplesPerm[,,1:33], cov90y[,,1:33], breaks=breaks)$stats[2,]
y95 = stats.bin(uSamplesPerm[,,1:33], cov95y[,,1:33], breaks=breaks)$stats[2,]
pdf("figures/coverageExperiment/cov50uDataDomain.pdf", width=5, height=5)
plot(x, y50, type="l", main="50% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.5, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov60uDataDomain.pdf", width=5, height=5)
plot(x, y60, type="l", main="60% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.6, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov70uDataDomain.pdf", width=5, height=5)
plot(x, y70, type="l", main="70% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.7, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov80uDataDomain.pdf", width=5, height=5)
plot(x, y80, type="l", main="80% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.8, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov90uDataDomain.pdf", width=5, height=5)
plot(x, y90, type="l", main="90% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.9, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov95uDataDomain.pdf", width=5, height=5)
plot(x, y95, type="l", main="95% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.95, lty=2)
dev.off()

## plot coverage versus u(s*) filtered to be outside of the data domain
# make the bins over values of u(s*)
# breaks = c(-Inf, seq(min(uSamples), max(uSamples), l=51)[-c(1, 50)], Inf)
breaks = qnorm(seq(1e-06, 1-1e-06, l=50))

# uSamples = array(dim=c(ns, n + m, nu))
# cov50 = array(dim=c(ns, nu, ny, m))
# uHats = array(dim=c(ns, m, nu, ny))
uSamplesPerm = aperm(uSamples, c(1, 3, 2))[,,-(1:n)]
out = stats.bin(uSamplesPerm[,,-(1:33)], cov50y[,,-(1:33)], breaks=breaks)
x = out$centers
y50 = out$stats[2,]
y60 = stats.bin(uSamplesPerm[,,-(1:33)], cov60y[,,-(1:33)], breaks=breaks)$stats[2,]
y70 = stats.bin(uSamplesPerm[,,-(1:33)], cov70y[,,-(1:33)], breaks=breaks)$stats[2,]
y80 = stats.bin(uSamplesPerm[,,-(1:33)], cov80y[,,-(1:33)], breaks=breaks)$stats[2,]
y90 = stats.bin(uSamplesPerm[,,-(1:33)], cov90y[,,-(1:33)], breaks=breaks)$stats[2,]
y95 = stats.bin(uSamplesPerm[,,-(1:33)], cov95y[,,-(1:33)], breaks=breaks)$stats[2,]
pdf("figures/coverageExperiment/cov50uExtrapolate.pdf", width=5, height=5)
plot(x, y50, type="l", main="50% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.5, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov60uExtrapolate.pdf", width=5, height=5)
plot(x, y60, type="l", main="60% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.6, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov70uExtrapolate.pdf", width=5, height=5)
plot(x, y70, type="l", main="70% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.7, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov80uExtrapolate.pdf", width=5, height=5)
plot(x, y80, type="l", main="80% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.8, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov90uExtrapolate.pdf", width=5, height=5)
plot(x, y90, type="l", main="90% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.9, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov95uExtrapolate.pdf", width=5, height=5)
plot(x, y95, type="l", main="95% Coverage", xlab="u(s)", ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.95, lty=2)
dev.off()

## do the same thing except with uHat instead of u
breaks = qnorm(seq(1e-06, 1-1e-06, l=50))

# uSamples = array(dim=c(ns, n + m, nu))
# cov50 = array(dim=c(ns, nu, ny, m))
# uHats = array(dim=c(ns, m, nu, ny))
uHatsY = apply(uHats, 1:3, mean)
uHatsPerm = aperm(uHatsY, c(1, 3, 2))
out = stats.bin(uHatsPerm, cov50y, breaks=breaks)
x = out$centers
y50 = out$stats[2,]
y60 = stats.bin(uHatsPerm, cov60y, breaks=breaks)$stats[2,]
y70 = stats.bin(uHatsPerm, cov70y, breaks=breaks)$stats[2,]
y80 = stats.bin(uHatsPerm, cov80y, breaks=breaks)$stats[2,]
y90 = stats.bin(uHatsPerm, cov90y, breaks=breaks)$stats[2,]
y95 = stats.bin(uHatsPerm, cov95y, breaks=breaks)$stats[2,]
pdf("figures/coverageExperiment/cov50uHat.pdf", width=5, height=5)
plot(x, y50, type="l", main="50% Coverage", xlab=TeX("$\\hat{u}(s)$"), ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.5, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov60uHat.pdf", width=5, height=5)
plot(x, y60, type="l", main="60% Coverage", xlab=TeX("$\\hat{u}(s)$"), ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.6, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov70uHat.pdf", width=5, height=5)
plot(x, y70, type="l", main="70% Coverage", xlab=TeX("$\\hat{u}(s)$"), ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.7, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov80uHat.pdf", width=5, height=5)
plot(x, y80, type="l", main="80% Coverage", xlab=TeX("$\\hat{u}(s)$"), ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.8, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov90uHat.pdf", width=5, height=5)
plot(x, y90, type="l", main="90% Coverage", xlab=TeX("$\\hat{u}(s)$"), ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.9, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov95uHat.pdf", width=5, height=5)
plot(x, y95, type="l", main="95% Coverage", xlab=TeX("$\\hat{u}(s)$"), ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.95, lty=2)
dev.off()

## Lot coverage versus uHat, filtered for only values with s* outside of the data domain
breaks = qnorm(seq(1e-06, 1-1e-06, l=50))

# uSamples = array(dim=c(ns, n + m, nu))
# cov50 = array(dim=c(ns, nu, ny, m))
# uHats = array(dim=c(ns, m, nu, ny))
uHatsY = apply(uHats, 1:3, mean)
uHatsPerm = aperm(uHatsY, c(1, 3, 2))
domainIndices=1:33
out = stats.bin(uHatsPerm[,,-domainIndices], cov50y[,,-domainIndices], breaks=breaks)
x = out$centers
y50 = out$stats[2,]
y60 = stats.bin(uHatsPerm[,,-domainIndices], cov60y[,,-domainIndices], breaks=breaks)$stats[2,]
y70 = stats.bin(uHatsPerm[,,-domainIndices], cov70y[,,-domainIndices], breaks=breaks)$stats[2,]
y80 = stats.bin(uHatsPerm[,,-domainIndices], cov80y[,,-domainIndices], breaks=breaks)$stats[2,]
y90 = stats.bin(uHatsPerm[,,-domainIndices], cov90y[,,-domainIndices], breaks=breaks)$stats[2,]
y95 = stats.bin(uHatsPerm[,,-domainIndices], cov95y[,,-domainIndices], breaks=breaks)$stats[2,]
pdf("figures/coverageExperiment/cov50uHatExtrapolate.pdf", width=5, height=5)
plot(x, y50, type="l", main="50% Coverage", xlab=TeX("$\\hat{u}(s)$"), ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.5, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov60uHatExtrapolate.pdf", width=5, height=5)
plot(x, y60, type="l", main="60% Coverage", xlab=TeX("$\\hat{u}(s)$"), ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.6, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov70uHatExtrapolate.pdf", width=5, height=5)
plot(x, y70, type="l", main="70% Coverage", xlab=TeX("$\\hat{u}(s)$"), ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.7, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov80uHatExtrapolate.pdf", width=5, height=5)
plot(x, y80, type="l", main="80% Coverage", xlab=TeX("$\\hat{u}(s)$"), ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.8, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov90uHatExtrapolate.pdf", width=5, height=5)
plot(x, y90, type="l", main="90% Coverage", xlab=TeX("$\\hat{u}(s)$"), ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.9, lty=2)
dev.off()

pdf("figures/coverageExperiment/cov95uHatExtrapolate.pdf", width=5, height=5)
plot(x, y95, type="l", main="95% Coverage", xlab=TeX("$\\hat{u}(s)$"), ylab="Coverage Probability", ylim=c(0, 1), 
     col="blue")
abline(h=.95, lty=2)
dev.off()

## Fix a given sample of s, s*, and u. Plot coverage over the distribution of y versus s*
# uSamples = array(dim=c(ns, n + m, nu))
# cov50 = array(dim=c(ns, nu, ny, m))
# uHats = array(dim=c(ns, m, nu, ny))

si = 1
ui = 1

pdf("figures/coverageExperiment/cov50uFixed.pdf", width=6, height=10)
par(mfrow=c(2,1))
thisRange = range(c(uSamples[si, -(1:n), ui], ySamples[si, , ui, ]))
plot(ss, uSamples[si, -(1:n), ui], type="l", col="blue", main="u(s)", 
     xlab="s", ylab="u(s)", ylim=thisRange)
# points(rep(sSamples[si,], ny), ySamples[si, , ui, ], pch=".")
rug(sSamples[si,])
plot(ss, cov50y[si, ui,], type="l", col="blue", main="50% Coverage (Fixed u, Sample Locations)", 
     xlab="s", ylab="Coverage probability")
abline(h=.5, lty=2)
rug(sSamples[si,])
dev.off()

pdf("figures/coverageExperiment/cov60uFixed.pdf", width=6, height=10)
par(mfrow=c(2,1))
thisRange = range(c(uSamples[si, -(1:n), ui], ySamples[si, , ui, ]))
plot(ss, uSamples[si, -(1:n), ui], type="l", col="blue", main="u(s)", 
     xlab="s", ylab="u(s)", ylim=thisRange)
# points(rep(sSamples[si,], ny), ySamples[si, , ui, ], pch=".")
rug(sSamples[si,])
plot(ss, cov60y[si, ui,], type="l", col="blue", main="60% Coverage (Fixed u, Sample Locations)", 
     xlab="s", ylab="Coverage probability")
abline(h=.6, lty=2)
rug(sSamples[si,])
dev.off()

pdf("figures/coverageExperiment/cov70uFixed.pdf", width=6, height=10)
par(mfrow=c(2,1))
thisRange = range(c(uSamples[si, -(1:n), ui], ySamples[si, , ui, ]))
plot(ss, uSamples[si, -(1:n), ui], type="l", col="blue", main="u(s)", 
     xlab="s", ylab="u(s)", ylim=thisRange)
# points(rep(sSamples[si,], ny), ySamples[si, , ui, ], pch=".")
rug(sSamples[si,])
plot(ss, cov70y[si, ui,], type="l", col="blue", main="70% Coverage (Fixed u, Sample Locations)", 
     xlab="s", ylab="Coverage probability")
abline(h=.7, lty=2)
rug(sSamples[si,])
dev.off()

pdf("figures/coverageExperiment/cov80uFixed.pdf", width=6, height=10)
par(mfrow=c(2,1))
thisRange = range(c(uSamples[si, -(1:n), ui], ySamples[si, , ui, ]))
plot(ss, uSamples[si, -(1:n), ui], type="l", col="blue", main="u(s)", 
     xlab="s", ylab="u(s)", ylim=thisRange)
# points(rep(sSamples[si,], ny), ySamples[si, , ui, ], pch=".")
rug(sSamples[si,])
plot(ss, cov80y[si, ui,], type="l", col="blue", main="80% Coverage (Fixed u, Sample Locations)", 
     xlab="s", ylab="Coverage probability")
abline(h=.8, lty=2)
rug(sSamples[si,])
dev.off()

pdf("figures/coverageExperiment/cov90uFixed.pdf", width=6, height=10)
par(mfrow=c(2,1))
thisRange = range(c(uSamples[si, -(1:n), ui], ySamples[si, , ui, ]))
plot(ss, uSamples[si, -(1:n), ui], type="l", col="blue", main="u(s)", 
     xlab="s", ylab="u(s)", ylim=thisRange)
# points(rep(sSamples[si,], ny), ySamples[si, , ui, ], pch=".")
rug(sSamples[si,])
plot(ss, cov90y[si, ui,], type="l", col="blue", main="90% Coverage (Fixed u, Sample Locations)", 
     xlab="s", ylab="Coverage probability")
abline(h=.9, lty=2)
rug(sSamples[si,])
dev.off()

pdf("figures/coverageExperiment/cov95uFixed.pdf", width=6, height=10)
par(mfrow=c(2,1))
thisRange = range(c(uSamples[si, -(1:n), ui], ySamples[si, , ui, ]))
plot(ss, uSamples[si, -(1:n), ui], type="l", col="blue", main="u(s)", 
     xlab="s", ylab="u(s)", ylim=thisRange)
# points(rep(sSamples[si,], ny), ySamples[si, , ui, ], pch=".")
rug(sSamples[si,])
plot(ss, cov95y[si, ui,], type="l", col="blue", main="95% Coverage (Fixed u, Sample Locations)", 
     xlab="s", ylab="Coverage probability")
abline(h=.95, lty=2)
rug(sSamples[si,])
dev.off()

## now make the 2d barplot of coverage versus u and nearest neighbor distance
# nnd = matrix(nrow=ns, ncol=m)
# uSamples = array(dim=c(ns, n + m, nu))
# cov50 = array(dim=c(ns, nu, ny, m))
# uHats = array(dim=c(ns, m, nu, ny))
# uGroups = inla.group(uSamples[,-(1:n),])
# nndGroups = array(inla.group(nnd), dim=c(ns, m, nu))
# nndGroups = aperm(nndGroups, c(1, 3, 2))
nndPerm = array(nnd, dim=c(ns, m, nu))
nndPerm = aperm(nndPerm, c(1, 3, 2))
uSamplesPerm = aperm(uSamples, c(1, 3, 2))
# binAverages = aggregate(c(cov50y), by=data.frame(uGroups=c(uGroups), nndGroups=c(nndGroups)), mean)

averagingFun = function(x) {
  if(length(x) < 10) {
    return(NA)
  } else {
    mean(x)
  }
}

par(mfrow=c(1,1))
cols=makeRedGrayBlueDivergingColors(64, c(0,1), .5, rev=FALSE)
png("figures/coverageExperiment/cov50uDist.png", width=500, height=500)
quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov50y), ny=30, nx=150, col=cols, 
           ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab="u(s)", main="Coverage Probability (50% Significance)", 
           FUN=averagingFun, zlim=c(0,1))
dev.off()

cols=makeRedGrayBlueDivergingColors(64, c(0,1), .6, rev=FALSE)
png("figures/coverageExperiment/cov60uDist.png", width=500, height=500)
quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov60y), ny=30, nx=150, col=cols, 
           ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab="u(s)", main="Coverage Probability (60% Significance)", 
           FUN=averagingFun, zlim=c(0,1))
dev.off()

cols=makeRedGrayBlueDivergingColors(64, c(0,1), .7, rev=FALSE)
png("figures/coverageExperiment/cov70uDist.png", width=500, height=500)
quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov70y), ny=30, nx=150, col=cols, 
           ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab="u(s)", main="Coverage Probability (70% Significance)", 
           FUN=averagingFun, zlim=c(0,1))
dev.off()

cols=makeRedGrayBlueDivergingColors(64, c(0,1), .8, rev=FALSE)
png("figures/coverageExperiment/cov80uDist.png", width=500, height=500)
quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov80y), ny=30, nx=150, col=cols, 
           ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab="u(s)", main="Coverage Probability (80% Significance)", 
           FUN=averagingFun, zlim=c(0,1))
dev.off()

cols=makeRedGrayBlueDivergingColors(64, c(0,1), .9, rev=FALSE)
png("figures/coverageExperiment/cov90uDist.png", width=500, height=500)
quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov90y), ny=30, nx=150, col=cols, 
           ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab="u(s)", main="Coverage Probability (90% Significance)", 
           FUN=averagingFun, zlim=c(0,1))
dev.off()

cols=makeRedGrayBlueDivergingColors(64, c(0,1), .95, rev=FALSE)
png("figures/coverageExperiment/cov95uDist.png", width=500, height=500)
quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov95y), ny=30, nx=150, col=cols, 
           ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab="u(s)", main="Coverage Probability (95% Significance)", 
           FUN=averagingFun, zlim=c(0,1))
dev.off()

## now make the 2d barplot of coverage versus uHat and nearest neighbor distance
# nnd = matrix(nrow=ns, ncol=m)
# uSamples = array(dim=c(ns, n + m, nu))
# cov50 = array(dim=c(ns, nu, ny, m))
# uHats = array(dim=c(ns, m, nu, ny))
nndPerm = array(nnd, dim=c(ns, m, nu, ny))
nndPerm = aperm(nndPerm, c(1, 3, 4, 2))
uHatsPerm = aperm(uHats, c(1, 3, 4, 2))

# set.seed(123)
downSampleYI = sample(1:ny, 30, replace = FALSE) # maximum allowable so far: 30. minimum not allowable so far: 50
uHatsPermDownSample = c(uHatsPerm[,,downSampleYI,])
nndPermDownSample = c(nndPerm[,,downSampleYI,])

par(mfrow=c(1,1))
cols=makeRedGrayBlueDivergingColors(64, c(0,1), .5, rev=FALSE)
png("figures/coverageExperiment/cov50uHatDist.png", width=500, height=500)
quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov50[,,downSampleYI,]), ny=30, nx=150, col=cols, 
           ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab=TeX("$\\hat{u}(s)$"), main="Coverage Probability (50% Significance)", 
           FUN=averagingFun, zlim=c(0,1))
dev.off()

cols=makeRedGrayBlueDivergingColors(64, c(0,1), .6, rev=FALSE)
png("figures/coverageExperiment/cov60uHatDist.png", width=500, height=500)
quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov60[,,downSampleYI,]), ny=30, nx=150, col=cols, 
           ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab=TeX("$\\hat{u}(s)$"), main="Coverage Probability (60% Significance)", 
           FUN=averagingFun, zlim=c(0,1))
dev.off()

cols=makeRedGrayBlueDivergingColors(64, c(0,1), .7, rev=FALSE)
png("figures/coverageExperiment/cov70uHatDist.png", width=500, height=500)
quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov70[,,downSampleYI,]), ny=30, nx=150, col=cols, 
           ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab=TeX("$\\hat{u}(s)$"), main="Coverage Probability (70% Significance)", 
           FUN=averagingFun, zlim=c(0,1))
dev.off()

cols=makeRedGrayBlueDivergingColors(64, c(0,1), .8, rev=FALSE)
png("figures/coverageExperiment/cov80uHatDist.png", width=500, height=500)
quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov80[,,downSampleYI,]), ny=30, nx=150, col=cols, 
           ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab=TeX("$\\hat{u}(s)$"), main="Coverage Probability (80% Significance)", 
           FUN=averagingFun, zlim=c(0,1))
dev.off()

cols=makeRedGrayBlueDivergingColors(64, c(0,1), .9, rev=FALSE)
png("figures/coverageExperiment/cov90uHatDist.png", width=500, height=500)
quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov90[,,downSampleYI,]), ny=30, nx=150, col=cols, 
           ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab=TeX("$\\hat{u}(s)$"), main="Coverage Probability (90% Significance)", 
           FUN=averagingFun, zlim=c(0,1))
dev.off()

cols=makeRedGrayBlueDivergingColors(64, c(0,1), .95, rev=FALSE)
png("figures/coverageExperiment/cov95uHatDist.png", width=500, height=500)
quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov95[,,downSampleYI,]), ny=30, nx=150, col=cols, 
           ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab=TeX("$\\hat{u}(s)$"), main="Coverage Probability (95% Significance)", 
           FUN=averagingFun, zlim=c(0,1))
dev.off()











