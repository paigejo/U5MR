# Code for INLA SPDE model:

# Make the mesh grid for INLA's SPDE spatial model.  Accounts for spherical coordinate system.  
# Note that max.n and globe don't seem to affect the grid, and by default 3779 vertices are 
# generated.
# doPlot: If doPlot is true, plots the generated mesh
# offset: argument to inla.mesh.create
# doGlobe: if TRUE, makes grid on the sphere, assumes coordinates are lat/lon.  
#          Otherwise, assumes coordinates are esting/northing.  If FALSE, a good default 
#          for locs is locs=projKenya(mort$lon, mort$lat)
# cutoff: minimum leg size of trianglation
getSPDEMeshGrid = function(locs=NULL, max.n=100, doPlot=TRUE, 
                           doGlobe=FALSE, offset=-.18, cutoff=ifelse(doGlobe, .15, 15)) {
  # base mesh off of actual DHS dataset unless user specifies otherwise
  if(is.null(locs)) {
    locs = cbind(mort$lon, mort$lat)
    if(!doGlobe)
      locs = projKenya(locs)
  }
  
  if(doGlobe) {
    # generate the mesh on the globe
    mesh = inla.mesh.create(locs, globe=10, refine=list(max.n=max.n), extend=list(offset=offset), 
                            cutoff=cutoff)
  }
  else {
    # generate mesh on R2
    mesh = inla.mesh.create(loc=locs, refine=list(max.n=max.n), extend=list(offset=offset), 
                            cutoff=cutoff)
  }
  
  # plot the mesh if user wants
  if(doPlot) {
    plot(mesh)
    plotMapDat(project=!doGlobe, border="blue")
  }
  
  mesh
}

# generate default priors for SPDE model
# from Lindgren Rue (2015) "Bayesian Spatial Modelling with R-INLA"
# sigma0: field standard deviation. Doesn't matter too much since it's the prior
getSPDEPrior = function(mesh, sigma0=1) {
  size <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
  range0 <- size/5
  kappa0 <- sqrt(8)/range0
  tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
  spde <- inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, +1),
                            B.kappa = cbind(log(kappa0), 0, -1), theta.prior.mean = c(0, 0),
                            theta.prior.prec = c(0.1, 1))
  
  spde
}

# fit the SPDE model given some data and some prediction locations
fitSPDEModel = function(obsCoords, obsNs=rep(25, nrow(obsCoords)), obsCounts, obsUrban, predCoords, predNs = rep(25, nrow(predCoords)), 
                        predUrban, prior=NULL, mesh=NULL, int.strategy="eb", strategy="gaussian", 
                        genCountyLevel=FALSE, popGrid=NULL, nPostSamples=100, kmRes=5, counties=as.character(unique(mort$admin1)), 
                        includeClustEffect=TRUE, verbose=TRUE, genRegionLevel=FALSE, regions=as.character(unique(mort$region)), 
                        keepPixelPreds=FALSE, genEALevel=FALSE, eaIndices=1:nrow(kenyaEAs), 
                        urbanEffect=TRUE) {
  
  # for enumeration area level predictions, set enumeration area end pixel level prediction indices
  eaMat = NULL
  if(genEALevel) {
    pixelIndices = (1:nrow(predCoords))[-eaIndices]
  }
  else {
    pixelIndices = 1:nrow(predCoords)
  }
  
  # make default mesh
  if(is.null(mesh)) {
    mesh = getSPDEMeshGrid(doPlot = FALSE)
  }
  
  # make default prior
  if(is.null(prior)) {
    prior = getSPDEPrior(mesh)
  }
  
  # make covariate matrices (intercept plus urban/rural depending on whether urban effect included)
  if(urbanEffect) {
    X = cbind(1, obsUrban)
    XPred = cbind(1, predUrban)
  }
  else {
    X = matrix(1, ncol=1, nrow=nrow(obsCoords))
    XPred = matrix(1, ncol=1, nrow=nrow(predCoords))
  }
  
  
  # construct A matrix for observations
  m = nrow(obsCoords)
  AEst = inla.spde.make.A(mesh, loc = obsCoords)
  
  # construct A matrix for predictions
  APred = inla.spde.make.A(mesh, loc = predCoords)
  
  # make inla stack
  ys = obsCounts
  n = ncol(AEst) # number of basis elements
  nObs = length(ys) # number of observations
  nPreds = nrow(predCoords)
  latticeInds = 1:n
  # clustIObs = (n+1):(n+nObs)
  # clustIPreds = (n+1):(n+nPreds)
  # clustIPreds = (n+nObs+1):(n+nObs+nPreds)
  # stack.est = inla.stack(A =list(AEst, 1, 1), 
  #                        effects =list(field=latticeInds, clust=clustIObs, X=X), 
  #                        data =list(y=ys, Ntrials = obsNs, link=1), 
  #                        tag ="est", 
  #                        remove.unused=FALSE)
  # stack.pred = inla.stack(A =list(APred, 1, 1),
  #                         effects =list(field=latticeInds, clust=clustIPreds, X=XPred),
  #                         data =list(y=NA, Ntrials = predNs, link=1),
  #                         tag ="pred",
  #                         remove.unused=FALSE)
  stack.est = inla.stack(A =list(AEst, 1), 
                         effects =list(field=latticeInds, X=X), 
                         data =list(y=ys, Ntrials = obsNs, link=1), 
                         tag ="est", 
                         remove.unused=FALSE)
  stack.pred = inla.stack(A =list(APred, 1),
                          effects =list(field=latticeInds, X=XPred),
                          data =list(y=NA, Ntrials = predNs, link=1),
                          tag ="pred",
                          remove.unused=FALSE)
  
  # make mesh index
  mesh.index <- inla.spde.make.index(name = "field", n.spde = prior$n.spde)
  
  # fit model
  control.inla = list(strategy=strategy, int.strategy=int.strategy) 
  allNs = c(obsNs, predNs)
  stack.full = inla.stack(stack.est, stack.pred, remove.unused=FALSE)
  stackDat = inla.stack.data(stack.full, spde=prior)
  # mod = inla(y ~ - 1 + X + f(field, model=prior) + f(clust, model="iid"), 
  #            data = stackDat, Ntrials=stackDat$Ntrials, 
  #            control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link), 
  #            family="binomial", verbose=TRUE, control.inla=control.inla, 
  #            control.compute=list(config=TRUE))
  mod = inla(y ~ - 1 + X + f(field, model=prior), 
             data = stackDat, Ntrials=stackDat$Ntrials, 
             control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link), 
             family="binomial", verbose=verbose, control.inla=control.inla, 
             control.compute=list(config=TRUE))
  
  # get predictive surface, SD, and data
  n = nrow(obsCoords)
  # obsInds = 1:n
  obsInds = inla.stack.index(stack.full, "est")$data
  predInds = inla.stack.index(stack.full, "pred")$data
  index = inla.stack.index(stack.full, "pred")$data
  linpred.mean = mod[["summary.fitted.values"]]$mean
  linpred.sd = mod[["summary.fitted.values"]]$sd
  
  # predInds = (n+1):(n+nrow(predCoords))
  preds = linpred.mean
  predSDs = linpred.sd
  
  countyPredMat = NULL
  if(genCountyLevel) {
    ## generate county level predictions
    
    # if not supplied, get grid of population densities for pop-weighted integration
    if(is.null(popGrid))
      popGrid = makeInterpPopGrid(kmRes=kmRes)
    
    # make sure prediction coordinates correspond to population density grid coordinates
    if(length(pixelIndices) != nrow(popGrid)) {
      stop("for county level predictions, prediction points must be on same grid as popGrid")
    }
    
    # generate samples from posterior
    postSamples = inla.posterior.sample(nPostSamples, mod)
    latentMat = sapply(postSamples, function(x) x$latent)
    
    # get logit mortality rates at the prediction locations
    index.pred <- inla.stack.index(stack.full,tag="pred")$data
    if(genEALevel) {
      # separate enumeration area and pixel level predictions
      indicesEA = index.pred[]
    }
    predMat <- latentMat[index.pred,]
    
    if(genEALevel) {
      # make sure to separate enumeration area and pixel level predictions
      eaMat = predMat[eaIndices, ]
      predMat = predMat[pixelIndices, ]
    }
    
    # the below code is commented since we removed the cluster effect from the model. The commented code is not yet tested
    # Also, if the cluster effect is included, it should only be included 4 the enumeration area predictions, 
    # not the areal level predictions
    # # remove cluster effect for integration
    # nPreds = nrow(predCoords)
    # clusterIndices = inla.stack.index(stack.full,tag="pred")$effects[(1 + mesh$n):(mesh$n + nPreds)]
    # clusterEffects = latentMat[clusterIndices, ]
    # noClustMat = predMat - clusterEffects
    
    # for each county, integrate predictions with respect to population density
    predCounties = popGrid$admin1
    pops = popGrid$popOrig
    integratePredsByCounty = function(countyName) {
      # subset data by county of interest
      countyI = predCounties == countyName
      countyCoords = predCoords[countyI,]
      countyPredMat = expit(predMat[countyI,])
      countyPops = pops[countyI]
      
      # compute average integral of predictions
      countyPops = countyPops/sum(countyPops)
      integrals = t(countyPredMat) %*% countyPops
    }
    countyPredMat <- t(sapply(counties, integratePredsByCounty))
  }
  
  regionPredMat = NULL
  if(genRegionLevel) {
    ## generate region level predictions
    
    # if not supplied, get grid of population densities for pop-weighted integration
    if(is.null(popGrid))
      popGrid = makeInterpPopGrid(kmRes=kmRes)
    
    # make sure prediction coordinates correspond to population density grid coordinates
    if(length(pixelIndices) != nrow(popGrid)) {
      stop("for county level predictions, prediction points must be on same grid as popGrid")
    }
    
    # generate samples from posterior
    if(is.null(predMat)) {
      postSamples = inla.posterior.sample(nPostSamples, mod)
      latentMat = sapply(postSamples, function(x) x$latent)
      
      # get logit mortality rates at the prediction locations
      index.pred <- inla.stack.index(stack.full,tag="pred")$data
      predMat <- latentMat[index.pred,]
      
      if(genEALevel) {
        # make sure to separate enumeration area and pixel level predictions
        eaMat = predMat[eaIndices, ]
        predMat = predMat[pixelIndices, ]
      }
    }
    
    # the below code is commented since we removed the cluster effect from the model. The commented code is not yet tested
    # # remove cluster effect for integration
    # nPreds = nrow(predCoords)
    # clusterIndices = inla.stack.index(stack.full,tag="pred")$effects[(1 + mesh$n):(mesh$n + nPreds)]
    # clusterEffects = latentMat[clusterIndices, ]
    # noClustMat = predMat - clusterEffects
    
    # for each county, integrate predictions with respect to population density
    predCounties = popGrid$admin1
    predRegions = countyToRegion(predCounties)
    pops = popGrid$popOrig
    integratePredsByRegion = function(regionName) {
      # subset data by county of interest
      regionI = predRegions == regionName
      regionCoords = predCoords[regionI,]
      regionPredMat = expit(predMat[regionI,])
      regionPops = pops[regionI]
      
      # compute average integral of predictions
      regionPops = regionPops/sum(regionPops)
      integrals = t(regionPredMat) %*% regionPops
    }
    regionPredMat <- t(sapply(regions, integratePredsByRegion))
  }
  
  # generate pixel level predictions if necessary
  pixelPreds = NULL
  if(keepPixelPreds) {
    pixelPreds = predMat
  }
  
  list(mod=mod, preds=preds, SDs=predSDs, obsInds=obsInds, predInds=index, mesh=mesh, 
       prior=prior, stack=stack.full, countyPredMat=countyPredMat, regionPredMat=regionPredMat, 
       pixelPredMat=pixelPreds, eaPredMat=eaMat)
}

testFitSPDEModel = function(predCoords=NULL, nPredPts=NULL, predUrban=NULL, seed=1234, numClusters=1000) {
  # generate simulated dataset
  sigmasq = .5^2
  effRange = 100
  beta0 = -1.5
  gamma = -.75
  tausq = .25^2
  out = simDat(kenyaEAs, beta0=beta0, margVar=sigmasq, tausq=tausq, gamma=gamma, effRange=effRange, 
               seed=seed, NC=40, nLayer=2, numClusters=numClusters)
  eaDat = out$eaDat
  clustDat = out$clustDat
  
  # load in simulated dataset
  # load("simDat.RData") # "eaDat"    "clustDat" "eaID"
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = rep(25, nrow(obsCoords))
  obsCounts = clustDat$died
  obsUrban = clustDat$urban
  
  # set default prediction locations and urbanicity
  if(is.null(predCoords)) {
    popGrid = makeInterpPopGrid()
    predCoords = cbind(popGrid$east, popGrid$north)
    nPredPts = nrow(predCoords)
    predUrban = popGrid$urban
  }
  
  sampleIs = sample(1:nrow(predCoords), nPredPts, FALSE)
  predCoords = predCoords[sampleIs,]
  predNs = rep(25, nPredPts)
  predUrban = predUrban[sampleIs]
  
  # fit model
  out = fitSPDEModel(obsCoords, obsNs=rep(25, nrow(obsCoords)), obsCounts, obsUrban, predCoords, predNs = rep(25, nrow(predCoords)), 
                     predUrban, genCountyLevel=TRUE, popGrid=popGrid, nPostSamples=100)
  mod = out$mod
  preds = out$preds
  sds = out$SDs
  obsInds = out$obsInds
  predInds = out$predInds
  mesh = out$mesh
  prior = out$prior
  stack = out$stack
  countyPredMat = out$countyPredMat
  
  # plot things
  # true datasets for all enerumeration areas and the cluster sample
  pdf("figures/spdeDataset.pdf", width=8, height=8)
  par(mfrow =c(2, 2))
  zlim = c(0, quantile(c(eaDat$died/eaDat$numChildren, clustDat$died/clustDat$numChildren, 
                         eaDat$trueProbDeath), probs=.975))
  quilt.plot(eaDat$east, eaDat$north, eaDat$died/eaDat$numChildren, main="All Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, obsCounts/obsNs, main="Sample Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(eaDat$east, eaDat$north, eaDat$trueProbDeath, main="All True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, clustDat$trueProbDeath, main="Sample True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  dev.off()
  
  # predictions
  pdf("figures/spdeTestPredictions.pdf", width=8, height=8)
  par(mfrow=c(2,2))
  zlim = range(c(0, preds[obsInds]))
  quilt.plot(obsCoords, obsCounts/obsNs, main="Observed Mortality", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, preds[obsInds], main="Posterior mean at observation locations", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(predCoords, preds[predInds], main="Posterior Mean", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(predCoords, sds[predInds], main="Prediction SD", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim)
  plotMapDat(project=TRUE)
  dev.off()
  
  # get difference between urban and rural mortality
  urbanDiff = mean(preds[obsInds][!obsUrban]) - mean(preds[obsInds][obsUrban])
  print(paste0("mean urban prediction difference: ", urbanDiff))
  
  # plot SPDE parameter posteriors
  result <- inla.spde2.result(mod, "field", prior)
  pdf("figures/spdeTestPosterior1.pdf", width=9, height=5)
  par(mfrow=c(1,2))
  plot(result[["marginals.range.nominal"]][[1]], type = "l",
       main = "Nominal range, posterior density")
  abline(v=effRange, col="green")
  
  plot(result$marginals.variance.nominal[[1]], type = "l",
       main = "Variance, posterior density")
  abline(v=sigmasq, col="green")
  dev.off()
  
  pdf("figures/spdeTestPosterior2.pdf", width=8, height=8)
  par(mfrow=c(2,2))
  plot(mod$marginals.fixed[[1]], type="l", main="Intercept posterior (logit scale)")
  abline(v=beta0, col="green")
  plot(mod$marginals.fixed[[2]], type="l", main="Urban effect posterior (logit scale)")
  abline(v=gamma, col="green")
  # plot(mod$, type="l", main="Urban effect posterior (logit scale)")
  # abline(v=tausq, col="green")
  clustMarg = inla.tmarginal(function(x) {1/exp(x)}, mod$internal.marginals.hyperpar$`Log precision for clust`)
  plot(clustMarg, type="l", main="Cluster variance posterior (logit scale)")
  abline(v=tausq, col="green")
  dev.off()
  
  countyMeans = rowMeans(countyPredMat)
  pdf("figures/spdeTestCountyPreds.pdf", width=5, height=5)
  plotMapDat(plotVar=countyMeans, adm1, project=TRUE, new=TRUE, main="Predicted mortality rate")
  dev.off()
  
  # print summary
  print(summary(mod))
  
  out
}

# same as testFitSPDEModel, accept generates data using spde simulation instead of
# LatticeKrig simulation
testFitSPDEModel2 = function(predCoords=NULL, nPredPts=NULL, predUrban=NULL, seed=1234, 
                             numClusters=1000, popGrid=NULL) {
  # generate simulated dataset
  sigmasq = .5^2
  effRange = 100
  beta0 = -1.5
  gamma = -.75
  tausq = .25^2
  out = simDat(kenyaEAs, beta0=beta0, margVar=sigmasq, tausq=tausq, gamma=gamma, effRange=effRange, 
               seed=seed, NC=40, nLayer=2, numClusters=numClusters)
  eaDat = out$eaDat
  clustDat = out$clustDat
  
  # load in simulated dataset
  # load("simDat.RData") # "eaDat"    "clustDat" "eaID"
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = rep(25, nrow(obsCoords))
  obsCounts = clustDat$died
  obsUrban = clustDat$urban
  
  # set default prediction locations and urbanicity
  if(is.null(predCoords)) {
    if(is.null(popGrid)) {
      popGrid = makeInterpPopGrid()
    }
    predCoords = cbind(popGrid$east, popGrid$north)
    nPredPts = nrow(predCoords)
    predUrban = popGrid$urban
  }
  
  sampleIs = sample(1:nrow(predCoords), nPredPts, FALSE)
  predCoords = predCoords[sampleIs,]
  predNs = rep(25, nPredPts)
  predUrban = predUrban[sampleIs]
  
  # fit model
  out = fitSPDEModel(obsCoords, obsNs=rep(25, nrow(obsCoords)), obsCounts, obsUrban, predCoords, predNs = rep(25, nrow(predCoords)), 
                     predUrban, genCountyLevel=TRUE, popGrid=popGrid, nPostSamples=100)
  mod = out$mod
  preds = out$preds
  sds = out$SDs
  obsInds = out$obsInds
  predInds = out$predInds
  mesh = out$mesh
  prior = out$prior
  stack = out$stack
  countyPredMat = out$countyPredMat
  
  # plot things
  # true datasets for all enerumeration areas and the cluster sample
  pdf("figures/spdeDataset.pdf", width=8, height=8)
  par(mfrow =c(2, 2))
  zlim = c(0, quantile(c(eaDat$died/eaDat$numChildren, clustDat$died/clustDat$numChildren, 
                         eaDat$trueProbDeath), probs=.975))
  quilt.plot(eaDat$east, eaDat$north, eaDat$died/eaDat$numChildren, main="All Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, obsCounts/obsNs, main="Sample Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(eaDat$east, eaDat$north, eaDat$trueProbDeath, main="All True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, clustDat$trueProbDeath, main="Sample True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  dev.off()
  
  # predictions
  pdf("figures/spdeTestPredictions.pdf", width=8, height=8)
  par(mfrow=c(2,2))
  zlim = range(c(0, preds[obsInds]))
  quilt.plot(obsCoords, obsCounts/obsNs, main="Observed Mortality", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, preds[obsInds], main="Posterior mean at observation locations", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(predCoords, preds[predInds], main="Posterior Mean", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(predCoords, sds[predInds], main="Prediction SD", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim)
  plotMapDat(project=TRUE)
  dev.off()
  
  # get difference between urban and rural mortality
  urbanDiff = mean(preds[obsInds][!obsUrban]) - mean(preds[obsInds][obsUrban])
  print(paste0("mean urban prediction difference: ", urbanDiff))
  
  # plot SPDE parameter posteriors
  result <- inla.spde2.result(mod, "field", prior)
  pdf("figures/spdeTestPosterior1.pdf", width=9, height=5)
  par(mfrow=c(1,2))
  plot(result[["marginals.range.nominal"]][[1]], type = "l",
       main = "Nominal range, posterior density")
  abline(v=effRange, col="green")
  
  plot(result$marginals.variance.nominal[[1]], type = "l",
       main = "Variance, posterior density")
  abline(v=sigmasq, col="green")
  dev.off()
  
  pdf("figures/spdeTestPosterior2.pdf", width=8, height=8)
  par(mfrow=c(2,2))
  plot(mod$marginals.fixed[[1]], type="l", main="Intercept posterior (logit scale)")
  abline(v=beta0, col="green")
  plot(mod$marginals.fixed[[2]], type="l", main="Urban effect posterior (logit scale)")
  abline(v=gamma, col="green")
  # plot(mod$, type="l", main="Urban effect posterior (logit scale)")
  # abline(v=tausq, col="green")
  clustMarg = inla.tmarginal(function(x) {1/exp(x)}, mod$internal.marginals.hyperpar$`Log precision for clust`)
  plot(clustMarg, type="l", main="Cluster variance posterior (logit scale)")
  abline(v=tausq, col="green")
  dev.off()
  
  countyMeans = rowMeans(countyPredMat)
  pdf("figures/spdeTestCountyPreds.pdf", width=5, height=5)
  plotMapDat(plotVar=countyMeans, adm1, project=TRUE, new=TRUE, main="Predicted mortality rate")
  dev.off()
  
  # print summary
  print(summary(mod))
  
  out
}

simSPDE = function(coords, nsim=1, mesh=NULL, effRange=(max(coords[,1])-min(coords[,1]))/3, margVar=1) {
  # generate mesh grid if necessary
  if(is.null(mesh)) {
    mesh = getSPDEMeshGrid(coords, doPlot = FALSE)
  }
  
  # calculate SPDE model parameters based on Lindgren Rue (2015) "Bayesian Spatial Modelling with R-INLA"
  meshSize <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
  # it is easier to use theta and set sigma0 to 1 then to set sigma0 and the effective range directly
  # kappa0 <- sqrt(8)/effRange * meshSize # since nu = 1
  # kappa0 <- sqrt(8)/effRange # since nu = 1
  # kappa0 = sqrt(8) / 5
  # logKappa = log(kappa0)
  sigma0 = 1
  # tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
  # logTau = log(tau0)
  
  # from page 5 of the paper listed above:
  logKappa = 0.5 * log(8)
  logTau = 0.5 * (lgamma(1) - (lgamma(2) + log(4*pi))) - logKappa
  theta = c(log(sqrt(margVar)), log(effRange))
  spde <- inla.spde2.matern(mesh, B.tau = cbind(logTau, -1, +1),
                            B.kappa = cbind(logKappa, 0, -1), theta.prior.mean = theta,
                            theta.prior.prec = c(0.1, 1))
  
  
  # generate A and Q precision matrix
  Q = inla.spde2.precision(spde, theta = theta)
  A = inla.spde.make.A(mesh, coords)
  
  # generate simulations
  simField = inla.qsample(nsim, Q)
  simDat = as.matrix(A %*% simField)
  
  simDat
}