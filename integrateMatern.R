integrateMatern = function(xlim=c(0,1), ylim=c(0,1), n=50, ranges=seq(.01, .1, l=100), nsim=100) {
  xs = seq(xlim[1], xlim[2], l=n)
  ys = seq(ylim[1], ylim[2], l=n)
  pts = make.surface.grid(list(x=xs, y=ys))
  integrals = numeric(length(ranges))
  for(i in 1:length(ranges)) {
    print(paste0("i: ", i))
    range = ranges[i]
    cov = stationary.cov(pts, theta=range)
    L = t(chol(cov))
    zsims = matrix(rnorm(nrow(pts)*nsim), nrow=nrow(pts), ncol=nsim)
    sims = L %*% zsims
    integrals[i] = mean(colSums(sims))
  }
  integrals
}

integrateMatern2 = function(xlim=c(0,1), ylim=c(0,1), n=50, ranges=seq(.01, .1, l=100)) {
  xs = seq(xlim[1], xlim[2], l=n)
  ys = seq(ylim[1], ylim[2], l=n)
  pts = make.surface.grid(list(x=xs, y=ys))
  integrals = numeric(length(ranges))
  for(i in 1:length(ranges)) {
    print(paste0("i: ", i))
    range = ranges[i]
    cov = stationary.cov(pts, theta=range)
    
    integrals[i] = mean(cov)
  }
  integrals
}