# this script plots results from the coverage experiment
library(fields)
library(RColorBrewer)
library(plotly)
source('~/git/U5MR/plotGenerator.R')

filled.contour2 <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...) 
  {
    # modification by Ian Taylor of the filled.contour function
    # to remove the key and facilitate overplotting with contour()
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
      stop("increasing 'x' and 'y' values expected")
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2]) * par("csi") * 2.54
    par(las = las)
    mar <- mar.orig
    plot.new()
    par(mar=mar)
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
      stop("no proper 'z' matrix specified")
    if (!is.double(z)) 
      storage.mode(z) <- "double"
    filled.contour(as.double(x), as.double(y), z, levels=as.double(levels), 
                            col = col)
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(x, side = 1)
        Axis(y, side = 2)
      }
    }
    else plot.axes
    if (frame.plot) 
      box()
    if (missing(plot.title)) 
      title(...)
    else plot.title
    invisible()
  }

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
out = stats.bin(uSamplesPerm[,,1:50], cov50y[,,1:50], breaks=breaks)
x = out$centers
y50 = out$stats[2,]
y60 = stats.bin(uSamplesPerm[,,1:50], cov60y[,,1:50], breaks=breaks)$stats[2,]
y70 = stats.bin(uSamplesPerm[,,1:50], cov70y[,,1:50], breaks=breaks)$stats[2,]
y80 = stats.bin(uSamplesPerm[,,1:50], cov80y[,,1:50], breaks=breaks)$stats[2,]
y90 = stats.bin(uSamplesPerm[,,1:50], cov90y[,,1:50], breaks=breaks)$stats[2,]
y95 = stats.bin(uSamplesPerm[,,1:50], cov95y[,,1:50], breaks=breaks)$stats[2,]
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
out = stats.bin(uSamplesPerm[,,-(1:50)], cov50y[,,-(1:50)], breaks=breaks)
x = out$centers
y50 = out$stats[2,]
y60 = stats.bin(uSamplesPerm[,,-(1:50)], cov60y[,,-(1:50)], breaks=breaks)$stats[2,]
y70 = stats.bin(uSamplesPerm[,,-(1:50)], cov70y[,,-(1:50)], breaks=breaks)$stats[2,]
y80 = stats.bin(uSamplesPerm[,,-(1:50)], cov80y[,,-(1:50)], breaks=breaks)$stats[2,]
y90 = stats.bin(uSamplesPerm[,,-(1:50)], cov90y[,,-(1:50)], breaks=breaks)$stats[2,]
y95 = stats.bin(uSamplesPerm[,,-(1:50)], cov95y[,,-(1:50)], breaks=breaks)$stats[2,]
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
domainIndices=1:50
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
    return(NA_real_)
  } else {
    mean(x)
  }
}

# use my.quilt.plot to average coverage values within each bin with at least ten observations
cols=makeRedBlueDivergingColors(64, c(0,1), .5, rev=FALSE)
test = my.quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov50y), ny=75, nx=150, FUN=averagingFun)
grid = test$grid
x = test$x
y = test$y
z = test$z
allPoints = expand.grid("u(s)"=x, Distance=y)
dataFrame = data.frame("u"=allPoints$`u(s)`, "NND"=allPoints$Distance, 
                       "Coverage"=c(z))

# now plot the results with contours
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(0, 1, by=.1))
png("figures/coverageExperiment/cov50uDist.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = u, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Coverage Probability (50% CI)") + xlab("u(s)") + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = u, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(0, 1, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$u), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()

# use my.quilt.plot to average coverage values within each bin with at least ten observations
cols=makeRedBlueDivergingColors(64, c(0,1), .6, rev=FALSE)
test = my.quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov60y), ny=75, nx=150, FUN=averagingFun, grid=grid)
x = test$x
y = test$y
z = test$z
allPoints = expand.grid("u(s)"=x, Distance=y)
dataFrame = data.frame("u"=allPoints$`u(s)`, "NND"=allPoints$Distance, 
                       "Coverage"=c(z))

# now plot the results with contours
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(0, 1, by=.1))
png("figures/coverageExperiment/cov60uDist.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = u, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Coverage Probability (60% CI)") + xlab("u(s)") + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = u, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(0, 1, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$u), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()

# use my.quilt.plot to average coverage values within each bin with at least ten observations
cols=makeRedBlueDivergingColors(64, c(0,1), .7, rev=FALSE)
test = my.quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov70y), ny=75, nx=150, FUN=averagingFun, grid=grid)
x = test$x
y = test$y
z = test$z
allPoints = expand.grid("u(s)"=x, Distance=y)
dataFrame = data.frame("u"=allPoints$`u(s)`, "NND"=allPoints$Distance, 
                       "Coverage"=c(z))

# now plot the results with contours
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(0, 1, by=.1))
png("figures/coverageExperiment/cov70uDist.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = u, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Coverage Probability (70% CI)") + xlab("u(s)") + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = u, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(0, 1, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$u), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()

# use my.quilt.plot to average coverage values within each bin with at least ten observations
cols=makeRedBlueDivergingColors(64, c(0,1), .8, rev=FALSE)
test = my.quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov80y), ny=75, nx=150, FUN=averagingFun, grid=grid)
x = test$x
y = test$y
z = test$z
allPoints = expand.grid("u(s)"=x, Distance=y)
dataFrame = data.frame("u"=allPoints$`u(s)`, "NND"=allPoints$Distance, 
                       "Coverage"=c(z))

# now plot the results with contours
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(0, 1, by=.1))
png("figures/coverageExperiment/cov80uDist.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = u, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Coverage Probability (80% CI)") + xlab("u(s)") + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = u, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(0, 1, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$u), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()

# use my.quilt.plot to average coverage values within each bin with at least ten observations
cols=makeRedBlueDivergingColors(64, c(0,1), .9, rev=FALSE)
test = my.quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov90y), ny=75, nx=150, FUN=averagingFun, grid=grid)
x = test$x
y = test$y
z = test$z
allPoints = expand.grid("u(s)"=x, Distance=y)
dataFrame = data.frame("u"=allPoints$`u(s)`, "NND"=allPoints$Distance, 
                       "Coverage"=c(z))

# now plot the results with contours
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(0, 1, by=.1))
png("figures/coverageExperiment/cov90uDist.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = u, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Coverage Probability (90% CI)") + xlab("u(s)") + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = u, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(0, 1, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$u), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()

# use my.quilt.plot to average coverage values within each bin with at least ten observations
cols=makeRedBlueDivergingColors(64, c(0,1), .95, rev=FALSE)
test = my.quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov95y), ny=75, nx=150, FUN=averagingFun, grid=grid)
x = test$x
y = test$y
z = test$z
allPoints = expand.grid("u(s)"=x, Distance=y)
dataFrame = data.frame("u"=allPoints$`u(s)`, "NND"=allPoints$Distance, 
                       "Coverage"=c(z))

# now plot the results with contours
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(0, 1, by=.1))
png("figures/coverageExperiment/cov95uDist.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = u, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Coverage Probability (95% CI)") + xlab("u(s)") + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = u, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(0, 1, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$u), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()

# par(mfrow=c(1,1))
# cols=makeRedGrayBlueDivergingColors(64, c(0,1), .5, rev=FALSE)
# png("figures/coverageExperiment/cov50uDist.png", width=500, height=500)
# my.quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov50y), ny=30, nx=150, col=cols, 
#            ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab="u(s)", main="Coverage Probability (50% Significance)", 
#            FUN=averagingFun, grid=grid, zlim=c(0,1))
# dev.off()
# 
# cols=makeRedGrayBlueDivergingColors(64, c(0,1), .6, rev=FALSE)
# png("figures/coverageExperiment/cov60uDist.png", width=500, height=500)
# my.quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov60y), ny=30, nx=150, col=cols, 
#            ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab="u(s)", main="Coverage Probability (60% Significance)", 
#            FUN=averagingFun, grid=grid, zlim=c(0,1))
# dev.off()
# 
# cols=makeRedGrayBlueDivergingColors(64, c(0,1), .7, rev=FALSE)
# png("figures/coverageExperiment/cov70uDist.png", width=500, height=500)
# my.quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov70y), ny=30, nx=150, col=cols, 
#            ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab="u(s)", main="Coverage Probability (70% Significance)", 
#            FUN=averagingFun, grid=grid, zlim=c(0,1))
# dev.off()
# 
# cols=makeRedGrayBlueDivergingColors(64, c(0,1), .8, rev=FALSE)
# png("figures/coverageExperiment/cov80uDist.png", width=500, height=500)
# my.quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov80y), ny=30, nx=150, col=cols, 
#            ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab="u(s)", main="Coverage Probability (80% Significance)", 
#            FUN=averagingFun, grid=grid, zlim=c(0,1))
# dev.off()
# 
# cols=makeRedGrayBlueDivergingColors(64, c(0,1), .9, rev=FALSE)
# png("figures/coverageExperiment/cov90uDist.png", width=500, height=500)
# my.quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov90y), ny=30, nx=150, col=cols, 
#            ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab="u(s)", main="Coverage Probability (90% Significance)", 
#            FUN=averagingFun, grid=grid, zlim=c(0,1))
# dev.off()
# 
# cols=makeRedGrayBlueDivergingColors(64, c(0,1), .95, rev=FALSE)
# png("figures/coverageExperiment/cov95uDist.png", width=500, height=500)
# my.quilt.plot(c(uSamplesPerm[,,-(1:n)]), c(nndPerm) / phi, c(cov95y), ny=30, nx=150, col=cols, 
#            ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab="u(s)", main="Coverage Probability (95% Significance)", 
#            FUN=averagingFun, grid=grid, zlim=c(0,1))
# dev.off()

## now make the 2d barplot of coverage versus uHat and nearest neighbor distance
# nnd = matrix(nrow=ns, ncol=m)
# uSamples = array(dim=c(ns, n + m, nu))
# cov50 = array(dim=c(ns, nu, ny, m))
# uHats = array(dim=c(ns, m, nu, ny))
nndPerm = array(nnd, dim=c(ns, m, nu, ny))
nndPerm = aperm(nndPerm, c(1, 3, 4, 2))
uHatsPerm = aperm(uHats, c(1, 3, 4, 2))

# set.seed(123)
# Sys.setenv("R_MAX_VSIZE"=8e9) # for avoiding the following error:
# Error: vector memory exhausted (limit reached?)
# 
# Enter a frame number, or 0 to exit   
# 
# 1: my.quilt.plot(uHatsPermDownSample, nndPermDownSample/phi, c(cov50[, , downSampleYI, ]), ny = 75, nx = 150, FUN = ave
#               2: as.image(z, x = x, nx = nx, ny = ny, grid = grid, FUN = FUN, na.rm = na.rm)
#               3: cbind(grid$x[temp$index[[1]]], grid$y[temp$index[[2]]])
#               4: cbind(grid$x[temp$index[[1]]], grid$y[temp$index[[2]]])
downSampleYI = sample(1:ny, 50, replace = FALSE) # maximum allowable so far: 50, but it takes ~5-10 minutes
downSampleYI = 1:100
uHatsPermDownSample = c(uHatsPerm[,,downSampleYI,])
nndPermDownSample = c(nndPerm[,,downSampleYI,])

# make sure there are at least 385 unique samples of u in each bin
uIs = cov50
count = 1
for(i in 1:dim(uIs)[1]) {
  for(j in 1:dim(uIs)[2]) {
    uIs[i,j,,] = count
    count = count + 1
  }
}
# .5/sqrt(385) * qnorm(.975)
# [1] 0.04994451

out <- my.quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(uIs[,,downSampleYI,]), ny=75, nx=150, 
                            getUniqueN = TRUE)
enoughUSamples = out$z >= 385
grid = list(x=out$x, y=out$y)

# use quilt.plot to average coverage values within each bin with
cols=makeRedBlueDivergingColors(64, c(0,1), .5, rev=FALSE)
# test = quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov50[,,downSampleYI,]), ny=75, nx=150, FUN=averagingFun, zlim=c(0,1))
test = my.quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov50[,,downSampleYI,]), ny=75, nx=150, grid=grid)
x = test$x
y = test$y
z = test$z
z[!enoughUSamples] = NA
# NOTE: could average overall samples by calling quilt.plot multiple times and averaging with respect to test$weights
allPoints = expand.grid(uHat=x, Distance=y)
dataFrame = data.frame("uHat"=allPoints$uHat, "NND"=allPoints$Distance, "Coverage"=c(z))

# now plot the results with contours
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(0.05, .95, by=.1))
png("figures/coverageExperiment/cov50uHatDist.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = uHat, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Coverage Probability (50% CI)") + xlab(TeX("$\\hat{u}(s)$")) + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = uHat, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(0.05, .95, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$uHat), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()

# use quilt.plot to average coverage values within each bin with at least ten observations
cols=makeRedBlueDivergingColors(64, c(0,1), .6, rev=FALSE)
test = my.quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov60[,,downSampleYI,]), ny=75, nx=150, grid=grid)
x = test$x
y = test$y
z = test$z
z[!enoughUSamples] = NA
# NOTE: could average overall samples by calling my.quilt.plot multiple times and averaging with respect to test$weights
allPoints = expand.grid(uHat=x, Distance=y)
dataFrame = data.frame("uHat"=allPoints$uHat, "NND"=allPoints$Distance, "Coverage"=c(z))

# now plot the results with contours
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(0.05, .95, by=.1))
png("figures/coverageExperiment/cov60uHatDist.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = uHat, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Coverage Probability (60% CI)") + xlab(TeX("$\\hat{u}(s)$")) + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = uHat, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(0.05, .95, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$uHat), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()

# use my.quilt.plot to average coverage values within each bin with at least ten observations
cols=makeRedBlueDivergingColors(64, c(0,1), .7, rev=FALSE)
test = my.quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov70[,,downSampleYI,]), ny=75, nx=150, grid=grid)
x = test$x
y = test$y
z = test$z
z[!enoughUSamples] = NA
# NOTE: could average overall samples by calling my.quilt.plot multiple times and averaging with respect to test$weights
allPoints = expand.grid(uHat=x, Distance=y)
dataFrame = data.frame("uHat"=allPoints$uHat, "NND"=allPoints$Distance, "Coverage"=c(z))

# now plot the results with contours
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(0.05, .95, by=.1))
png("figures/coverageExperiment/cov70uHatDist.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = uHat, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Coverage Probability (70% CI)") + xlab(TeX("$\\hat{u}(s)$")) + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = uHat, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(0.05, .95, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$uHat), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()

# use my.quilt.plot to average coverage values within each bin with at least ten observations
cols=makeRedBlueDivergingColors(64, c(0,1), .8, rev=FALSE)
test = my.quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov80[,,downSampleYI,]), ny=75, nx=150, grid=grid)
x = test$x
y = test$y
z = test$z
z[!enoughUSamples] = NA
# NOTE: could average overall samples by calling my.quilt.plot multiple times and averaging with respect to test$weights
allPoints = expand.grid(uHat=x, Distance=y)
dataFrame = data.frame("uHat"=allPoints$uHat, "NND"=allPoints$Distance, "Coverage"=c(z))

# now plot the results with contours
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(0.05, .95, by=.1))
png("figures/coverageExperiment/cov80uHatDist.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = uHat, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Coverage Probability (80% CI)") + xlab(TeX("$\\hat{u}(s)$")) + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = uHat, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(0.05, .95, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$uHat), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()

# use my.quilt.plot to average coverage values within each bin with at least ten observations
cols=makeRedBlueDivergingColors(64, c(0,1), .9, rev=FALSE)
test = my.quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov90[,,downSampleYI,]), ny=75, nx=150, grid=grid)
x = test$x
y = test$y
z = test$z
z[!enoughUSamples] = NA
# NOTE: could average overall samples by calling my.quilt.plot multiple times and averaging with respect to test$weights
allPoints = expand.grid(uHat=x, Distance=y)
dataFrame = data.frame("uHat"=allPoints$uHat, "NND"=allPoints$Distance, "Coverage"=c(z))

# now plot the results with contours
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(0.05, .95, by=.1))
png("figures/coverageExperiment/cov90uHatDist.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = uHat, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Coverage Probability (90% CI)") + xlab(TeX("$\\hat{u}(s)$")) + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = uHat, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(0.05, .95, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$uHat), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()

# use my.quilt.plot to average coverage values within each bin with at least ten observations
cols=makeRedBlueDivergingColors(64, c(0,1), .95, rev=FALSE)
test = my.quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov95[,,downSampleYI,]), ny=75, nx=150, grid=grid)
x = test$x
y = test$y
z = test$z
z[!enoughUSamples] = NA
# NOTE: could average overall samples by calling my.quilt.plot multiple times and averaging with respect to test$weights
allPoints = expand.grid(uHat=x, Distance=y)
dataFrame = data.frame("uHat"=allPoints$uHat, "NND"=allPoints$Distance, "Coverage"=c(z))

# now plot the results with contours
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(0, 1, by=.1))
png("figures/coverageExperiment/cov95uHatDist.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = uHat, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Coverage Probability (95% CI)") + xlab(TeX("$\\hat{u}(s)$")) + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = uHat, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(0, 1, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$uHat), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()

# par(mfrow=c(1,1))
# cols=makeRedGrayBlueDivergingColors(64, c(0,1), .5, rev=FALSE)
# png("figures/coverageExperiment/cov50uHatDist.png", width=500, height=500)
# my.quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov50[,,downSampleYI,]), ny=30, nx=150, col=cols, 
#            ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab=TeX("$\\hat{u}(s)$"), main="Coverage Probability (50% Significance)", 
#            FUN=averagingFun, grid=grid)
# dev.off()

# cols=makeRedGrayBlueDivergingColors(64, c(0,1), .6, rev=FALSE)
# png("figures/coverageExperiment/cov60uHatDist.png", width=500, height=500)
# my.quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov60[,,downSampleYI,]), ny=30, nx=150, col=cols, 
#            ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab=TeX("$\\hat{u}(s)$"), main="Coverage Probability (60% Significance)", 
#            FUN=averagingFun, grid=grid)
# dev.off()
# 
# cols=makeRedGrayBlueDivergingColors(64, c(0,1), .7, rev=FALSE)
# png("figures/coverageExperiment/cov70uHatDist.png", width=500, height=500)
# my.quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov70[,,downSampleYI,]), ny=30, nx=150, col=cols, 
#            ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab=TeX("$\\hat{u}(s)$"), main="Coverage Probability (70% Significance)", 
#            FUN=averagingFun, grid=grid)
# dev.off()
# 
# cols=makeRedGrayBlueDivergingColors(64, c(0,1), .8, rev=FALSE)
# png("figures/coverageExperiment/cov80uHatDist.png", width=500, height=500)
# my.quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov80[,,downSampleYI,]), ny=30, nx=150, col=cols, 
#            ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab=TeX("$\\hat{u}(s)$"), main="Coverage Probability (80% Significance)", 
#            FUN=averagingFun, grid=grid)
# dev.off()
# 
# cols=makeRedGrayBlueDivergingColors(64, c(0,1), .9, rev=FALSE)
# png("figures/coverageExperiment/cov90uHatDist.png", width=500, height=500)
# my.quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov90[,,downSampleYI,]), ny=30, nx=150, col=cols, 
#            ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab=TeX("$\\hat{u}(s)$"), main="Coverage Probability (90% Significance)", 
#            FUN=averagingFun, grid=grid)
# dev.off()
# 
# cols=makeRedGrayBlueDivergingColors(64, c(0,1), .95, rev=FALSE)
# png("figures/coverageExperiment/cov95uHatDist.png", width=500, height=500)
# my.quilt.plot(uHatsPermDownSample, nndPermDownSample / phi, c(cov95[,,downSampleYI,]), ny=30, nx=150, col=cols, 
#            ylab="(Nearest Neighbor Distance)/(Correlation Range)", xlab=TeX("$\\hat{u}(s)$"), main="Coverage Probability (95% Significance)", 
#            FUN=averagingFun, grid=grid)
# dev.off()

# get theoretical coverage (averaging over y)
# getCoverageFixedU = function(u0=0, us, s, sigma=1, tau=.1, kappa=1, alpha=.05) {
#   rho = stationary.cov(matrix(c(0, s), ncol=1), Covariance="Matern", theta=1, smoothness=kappa)[1,2]
#   center = (sigma^2 + tau^2)/(rho * sigma^2) * (us - u0) / tau
#   sd = ((1 + tau^2 / sigma^2) / rho)  *  (sigma / tau)  *  sqrt(1 - (rho^2 * sigma^2) / (sigma^2 + tau^2))
#   z = qnorm(1-alpha/2)
#   pnorm(center + z * sd) - pnorm(center - z * sd)
# }
# 
# # integrate getCoverage over initial value of u(0)
# getCoverageRandomUConditional = function(us, s, sigma=1, tau=.1, kappa=1, alpha=.05) {
#   # get conditional distribution of u(0) given u(s)
#   rho = stationary.cov(matrix(c(0, s), ncol=1), Covariance="Matern", theta=1, smoothness=kappa)[1,2]
#   center = rho * us
#   sd = sqrt(sigma^2 - rho^2 * sigma^2)
#   
#   integrand = function(u0) {
#     getCoverageFixedU(u0, us=us, s=s, sigma=sigma, tau=tau, kappa=kappa, alpha=alpha) * dnorm(center, sd=sd)
#   }
#   out = integrate(integrand, center - 5*sd, center + 5*sd)
#   out$value
# }
# 
# getCoverageRandomUMarginal = function(us, s, sigma=1, tau=.1, kappa=1, alpha=.05) {
#   # get marginal distribution of u(0)
#   center = 0
#   sd = sigma
#   
#   integrand = function(u0) {
#     getCoverageFixedU(u0, us=us, s=s, sigma=sigma, tau=tau, kappa=kappa, alpha=alpha) * dnorm(center, sd=sd)
#   }
#   out = integrate(integrand, center - 5*sd, center + 5*sd)
#   out$value
# }
# 
# 
# getCoverageU = function(us, s, sigma=1, tau=.1, kappa=1, alpha=.05) {
#   covMat = stationary.cov(matrix(c(0, s), ncol=1), Covariance="Matern", theta=1, smoothness=kappa)
#   rho = covMat[1,2:ncol(covMat)]
#   # center = (sigma^2 + tau^2)/(rho * sigma^2) * us / sqrt(tau^2 + sigma^2)
#   # sd = ((1 + tau^2 / sigma^2) / rho)  *  (sigma / sqrt(tau^2 + sigma^2))  *  sqrt(1 - (rho^2 * sigma^2) / (sigma^2 + tau^2))
#   center = sqrt(sigma^2 + tau^2)/(rho * sigma^2) * us
#   sd = sqrt(sigma^2 + tau^2) / (rho * sigma^2) * sigma * sqrt(1 - rho^2 * sigma^2 / (sigma^2 + tau^2))
#   z = qnorm(1-alpha/2)
#   pnorm(center + z * sd) - pnorm(center - z * sd)
# }

# get coverage over the distribution of y(0) | u(s)
getCoverageU = function(us, s, sigma=1, tau=.1, kappa=1, alpha=.05) {
  covMat = stationary.cov(matrix(c(0, s), ncol=1), Covariance="Matern", theta=1, smoothness=kappa)
  rho = covMat[1,2:ncol(covMat)]
  
  center = sqrt((1 - rho^2) * sigma^2 + tau^2) / (rho * sigma^2)  *  us
  sd = sqrt(sigma^2 + tau^2) / (rho * sigma)
  z = qnorm(1-alpha/2)
  pnorm(center + z * sd) - pnorm(center - z * sd)
}

uGrid = seq(-4, 4, l=150)
NNDGrid= seq(0, 10, l=100)
theoreticalCoverages95 = outer(uGrid, NNDGrid, FUN=getCoverageU)
theoreticalCoverages9 = outer(uGrid, NNDGrid, FUN=getCoverageU, alpha=.1)
theoreticalCoverages8 = outer(uGrid, NNDGrid, FUN=getCoverageU, alpha=.2)
theoreticalCoverages7 = outer(uGrid, NNDGrid, FUN=getCoverageU, alpha=.3)
theoreticalCoverages6 = outer(uGrid, NNDGrid, FUN=getCoverageU, alpha=.4)
theoreticalCoverages5 = outer(uGrid, NNDGrid, FUN=getCoverageU, alpha=.5)
allPoints = expand.grid(u=uGrid, Distance=NNDGrid)
dataFrame = data.frame("u"=allPoints$u, "NND"=allPoints$Distance, "Coverage"=c(theoreticalCoverages95))

# now plot the results with contours
cols=makeRedBlueDivergingColors(64, c(0,1), .95, rev=FALSE)
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(.1, .9, by=.1))
png("figures/coverageExperiment/cov95uDistTheoretical1point.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = u, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Theoretical Coverage Probability (95% CI)") + xlab("u(s)") + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = u, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(.1, .9, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$u), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()

allPoints = expand.grid(u=uGrid, Distance=NNDGrid)
dataFrame = data.frame("u"=allPoints$u, "NND"=allPoints$Distance, "Coverage"=c(theoreticalCoverages9))
cols=makeRedBlueDivergingColors(64, c(0,1), .9, rev=FALSE)
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(.1, .9, by=.1))
png("figures/coverageExperiment/cov90uDistTheoretical1point.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = u, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Theoretical Coverage Probability (90% CI)") + xlab("u(s)") + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = u, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(.1, .9, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$u), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()

allPoints = expand.grid(u=uGrid, Distance=NNDGrid)
dataFrame = data.frame("u"=allPoints$u, "NND"=allPoints$Distance, "Coverage"=c(theoreticalCoverages8))
cols=makeRedBlueDivergingColors(64, c(0,1), .8, rev=FALSE)
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(.1, .9, by=.1))
png("figures/coverageExperiment/cov80uDistTheoretical1point.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = u, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Theoretical Coverage Probability (80% CI)") + xlab("u(s)") + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = u, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(.1, .9, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$u), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()

allPoints = expand.grid(u=uGrid, Distance=NNDGrid)
dataFrame = data.frame("u"=allPoints$u, "NND"=allPoints$Distance, "Coverage"=c(theoreticalCoverages7))
cols=makeRedBlueDivergingColors(64, c(0,1), .7, rev=FALSE)
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(.1, .9, by=.1))
png("figures/coverageExperiment/cov70uDistTheoretical1point.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = u, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Theoretical Coverage Probability (70% CI)") + xlab("u(s)") + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = u, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(.1, .9, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$u), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()

allPoints = expand.grid(u=uGrid, Distance=NNDGrid)
dataFrame = data.frame("u"=allPoints$u, "NND"=allPoints$Distance, "Coverage"=c(theoreticalCoverages6))
cols=makeRedBlueDivergingColors(64, c(0,1), .6, rev=FALSE)
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(.1, .9, by=.1))
png("figures/coverageExperiment/cov60uDistTheoretical1point.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = u, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Theoretical Coverage Probability (60% CI)") + xlab("u(s)") + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = u, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(.1, .9, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$u), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()

allPoints = expand.grid(u=uGrid, Distance=NNDGrid)
dataFrame = data.frame("u"=allPoints$u, "NND"=allPoints$Distance, "Coverage"=c(theoreticalCoverages5))
cols=makeRedBlueDivergingColors(64, c(0,1), .5, rev=FALSE)
myPalette <- colorRampPalette(cols)
sc <- scale_fill_gradientn(colours = myPalette(64), limits=c(0, 1), breaks=seq(.1, .9, by=.1))
png("figures/coverageExperiment/cov50uDistTheoretical1point.png", width=500, height=500)
ggplot() + geom_raster(data = dataFrame, aes(x = u, y = NND, fill = Coverage), interpolate=FALSE) + 
  sc + ggtitle("Theoretical Coverage Probability (50% CI)") + xlab("u(s)") + ylab("(Nearest Neighbor Distance)/(Correlation Range)") + 
  stat_contour(data = dataFrame, aes(x = u, y = NND, z=Coverage), colour="black", size=.25, show.legend=TRUE, breaks=seq(.1, .9, by=.1)) + 
  scale_x_continuous(limits = range(dataFrame$u), expand=c(0,0)) + 
  scale_y_continuous(limits = range(dataFrame$NND), expand=c(0,0)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = guide_colourbar(barheight = unit( 3 , "in" ),
                                ticks.colour = "black",
                                ticks.linewidth = 1))
dev.off()





# out = runCompareModels2(maxDataSets = 2, saveResults=FALSE)