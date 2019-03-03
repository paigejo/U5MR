# Kenya utility scripts

# for plotting administration data
# project: if FALSE, plot with lon/lat coordinates.  Otherwise, plot with projected coords 
#          using projKenya function.  This can be used when plotting the projected `east' 
#          and `north' variables in kenyaEAs for instance.
# ...: arguments to polygon function
plotMapDat = function(mapDat = adm0, plotVar=NULL, varCounties=as.character(unique(mort$admin1)), zlim=NULL, project=FALSE, cols=tim.colors(), 
                      legend.mar=6, new=FALSE, plotArgs=NULL, main=NULL, xlim=NULL, xlab=NULL, 
                      ylim=NULL, ylab=NULL, ...) {
  # do setup for ploting data by county if necessary
  if(!is.null(plotVar)) {
    if(is.null(zlim)) {
      zlim = range(plotVar)
    }
    
    # get region names from map data
    regionNames = mapDat@data$NAME_1
    
    # make sure county names are consistent for mapDat == adm1
    regionNames[regionNames == "Elgeyo-Marakwet"] = "Elgeyo Marakwet"
    regionNames[regionNames == "Trans Nzoia"] = "Trans-Nzoia"
  }
  
  # generate new plot if necessary
  if(new) {
    # set graphical parameters so the legend won't overlap with plot
    currPar = par()
    newPar = currPar
    newMar = newPar$mar
    newMar[4] = max(newMar[4], legend.mar)
    newPar$mar = newMar
    if(currPar$mar[4] != newMar[4])
      suppressWarnings({par(newPar)})
    
    if(is.null(plotArgs)) {
      if(project) {
        if(is.null(xlab))
          xlab = "East (km)"
        if(is.null(xlim))
          xlim = eastLim
        if(is.null(ylab))
          ylab = "North (km)"
        if(is.null(ylim))
          ylim = northLim
      }
      else {
        if(is.null(xlab))
          xlab = "Longitude"
        if(is.null(xlim))
          xlim = kenyaLonRange
        if(is.null(ylab))
          ylab = "Latitude"
        if(is.null(ylim))
          ylim = kenyaLatRange
      }
      if(is.null(main))
        main = ""
      plotArgs = list(main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
    }
    do.call("plot", c(list(1, 2, type="n"), plotArgs))
  }
  
  # add polygons to plot
  polys = mapDat@polygons
  plotCounty = function(i) {
    countyPolys = polys[[i]]@Polygons
    
    if(is.null(plotVar)) {
      if(!project)
        sapply(1:length(countyPolys), function(x) {do.call("polygon", c(list(countyPolys[[x]]@coords), list(...)))})
      else
        sapply(1:length(countyPolys), function(x) {do.call("polygon", c(list(projKenya(countyPolys[[x]]@coords)), list(...)))})
    }
    else {
      # get index of plotVar corresponding to this county
      thisI = which(varCounties == regionNames[i])
      
      # get color to plot
      vals = c(zlim, plotVar[thisI])
      vals = vals-vals[1]
      vals = vals/(vals[2] - vals[1])
      col = cols[round(vals[3]*(length(cols)-1))+1]
      
      if(!project)
        sapply(1:length(countyPolys), function(x) {do.call("polygon", c(list(countyPolys[[x]]@coords, col=col), list(...)))})
      else
        sapply(1:length(countyPolys), function(x) {do.call("polygon", c(list(projKenya(countyPolys[[x]]@coords), col=col), list(...)))})
    }
    
  }
  sapply(1:length(polys), plotCounty)
  
  if(!is.null(plotVar)) {
    # add legend
    image.plot(zlim=zlim, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE, 
               col=cols, add = TRUE)
  }
  invisible(NULL)
}

# for computing what administrative regions the given points are in
# project: project to longitude/latitude coordinates
getRegion = function(points, mapDat = adm1, project=FALSE) {
  # project points to lon/lat coordinate system if user specifies
  if(project)
    points = projKenya(points, inverse=TRUE)
  
  regionNames = mapDat@data$NAME_1
  
  # make sure county names are consistent for mapDat == adm1
  regionNames[regionNames == "Elgeyo-Marakwet"] = "Elgeyo Marakwet"
  regionNames[regionNames == "Trans Nzoia"] = "Trans-Nzoia"
  
  # get region map polygons and set helper function for testing if points are in the regions
  polys = mapDat@polygons
  inRegion = function(i) {
    countyPolys = polys[[i]]@Polygons
    inside = sapply(1:length(countyPolys), function(x) {in.poly(points, countyPolys[[x]]@coords, inflation=0)})
    insideAny = apply(inside, 1, any)
    return(insideAny*i)
  }
  out = sapply(1:length(polys), inRegion)
  multipleRegs = apply(out, 1, function(vals) {sum(vals != 0) > 1})
  regionID = apply(out, 1, function(vals) {match(1, vals != 0)})
  regionNameVec = regionNames[regionID]
  list(regionID=regionID, regionNames=regionNameVec, multipleRegs=multipleRegs)
}

# convert county to region
countyToRegion = function(countyNames) {
  regionIs = match(as.character(countyNames), as.character(ctp[,1]))
  as.character(ctp[regionIs,2])
}

# do stratified sampling
# stratumVec: vector of stratum values, and each row being a sample
# sizes: vector with same length as unique(stratumVec) indicating the number of indices to sample
# probs: within each unique stratum level sample with probability proportional to probs (e.g. population density)
sampleByStratum = function(stratumVec, sizes, probs=1/length(stratumVec)) {
  uniqueStratumVec = unique(stratumVec)
  probMat = sapply(1:length(sizes), function(x) {probs*(stratumVec == uniqueStratumVec[x])/sum(probs*(stratumVec == uniqueStratumVec[x]))})
  unlist(sapply(1:length(sizes), function(x) {sample(1:length(stratumVec), size=sizes[x], prob=probMat[,x])}))
}

# do multilayer stratified sampling
# strataMat: each column is a stratum and each row is a sample
# sizes: vector of sizes of length equal to the number of unique 
#        combinations of values of levels in strataMat.  Number of 
#        elements to sample from that combination of strata levels
# probs: sample proportionally to probs within each strata level 
#        combination
sampleByStrata = function(strataMat, sizes, probs=1/nrow(strataMat)) {
  if(ncol(strataMat) == 1) 
    stop("only one stratum was input.  Use sampleByStratum in this case.")
  
  strataVec = strataMat[,1]
  for(i in 2:ncol(strataMat)) {
    strataVec = interaction(strataVec, strataMat[,i], drop=TRUE)
  }
  
  samples = sampleByStratum(strataVec, sizes, probs)
  
  list(strataVec=strataVec, samples=samples)
}

# generate a grid with a fixed latitude and longitude resolution.
# either set degrees per cell (res), number of grid cells on longest axis (nc), 
# or number of grid cells on both axes (nx and ny).
getKenyaGrid = function(res=.25, nc=NULL, nx=NULL, ny=NULL) {
  kenyaLonRange = c(33.5, 42)
  kenyaLatRange = c(-5,5.5)
  lonLen = kenyaLonRange[2] - kenyaLonRange[1]
  latLen = kenyaLatRange[2] - kenyaLatRange[1]
  
  # set individual axis grids
  if(!is.null(nc))
    res = latLen/nc
  if(is.null(nx) && is.null(ny)) {
    # we must go off of res
    lons = seq(kenyaLonRange[1], kenyaLonLength[2], by=res)
    lats = seq(kenyaLatRange[1], kenyaLatLength[2], by=res)
  }
  else {
    lons = seq(kenyaLonRange[1], kenyaLonLength[2], l=nx)
    lats = seq(kenyaLatRange[1], kenyaLatLength[2], l=ny)
  }
  
  # get full grid
  fullGrid = make.surface.grid(list(lon=lons, lat=lats))
  
  # return subset of points inside Kenya polygon
  polys = adm0@polygons
  kenyaPoly = polys[[1]]@Polygons[[77]]@coords
  inKenya = in.poly(fullGrid, kenyaPoly)
  fullGrid[inKenya,]
}

##### put Kenya population density on a grid of the chosen resolution
makeInterpPopGrid = function(kmRes=5) {
  # load population density data
  require(raster)
  
  # pop = raster("Kenya2014Pop/worldpop_total_1y_2014_00_00.tif", values= TRUE)
  load("Kenya2014Pop/pop.RData")
  
  # get a rectangular grid
  eastGrid = seq(eastLim[1], eastLim[2], by=kmRes)
  northGrid = seq(northLim[1], northLim[2], by=kmRes)
  utmGrid = make.surface.grid(list(east=eastGrid, north=northGrid))
  
  # project coordinates into lat/lon
  lonLatGrid = projKenya(utmGrid, inverse=TRUE)
  
  # subset grid so it's in Kenya
  polys = adm0@polygons
  kenyaPoly = polys[[1]]@Polygons[[77]]@coords
  inKenya = in.poly(lonLatGrid, kenyaPoly)
  utmGrid = utmGrid[inKenya,]
  lonLatGrid = lonLatGrid[inKenya,]
  
  # get population density at those coordinates
  interpPopVals = extract(pop, SpatialPoints(lonLatGrid),method="bilinear")
  
  # compute counties associated with locations
  counties = getRegion(lonLatGrid, adm1)$regionNames
  
  # determine which points are urban
  newPop = data.frame(list(lon=lonLatGrid[,1], lat=lonLatGrid[,2], popOrig=interpPopVals, admin1=counties))
  threshes = setThresholds2()
  popThreshes = sapply(1:nrow(newPop), function(i) {threshes$threshes[threshes$counties == newPop$admin1[i]]})
  urban = newPop$popOrig > popThreshes
  newPop$urban = urban
  
  newPop$east = utmGrid[,1]
  newPop$north = utmGrid[,2]
  
  newPop
}

##### some random convenient function for breaking down number of urban counties for simulated 
##### and real datasets
numPerCounty = function(datSet, byUrban=TRUE, counties=clustpc$County) {
  if(!byUrban)
    results = aggregate(datSet$urban, list(datSet$admin1), length)
  else
    results = aggregate(datSet$urban, list(datSet$admin1), sum)
  sortI = match(counties, results$Group.1)
  out = matrix(results$x[sortI], ncol=1)
  rownames(out) = results$Group.1[sortI]
  out
}

propUrbPerCounty = function(datSet, counties=clustpc$County) {
  numPerCounty(datSet)/numPerCounty(datSet, FALSE, counties)
}

##### project from lat/lon to UTM northing/easting in kilometers.  Use epsg=21097
# either pass lon/east and lat/north, or a matrix with 2 columns: first being lon/east, second being lat/north
# inverse: if FALSE, projects from lon/lat to easting/northing.  Else from easting/northing to lon/lat
projKenya = function(lon, lat=NULL, inverse=FALSE) {
  if(is.null(lat)) {
    lat = lon[,2]
    lon = lon[,1]
  }
  
  if(!inverse) {
    # from lon/lat coords to easting/northing
    lonLatCoords = SpatialPoints(cbind(lon, lat), proj4string=CRS("+proj=longlat"))
    coordsUTM = spTransform(lonLatCoords, CRS("+init=epsg:21097 +units=km"))
    out = attr(coordsUTM, "coords")
  }
  else {
    # from easting/northing coords to lon/lat
    east = lon
    north = lat
    coordsUTM = SpatialPoints(cbind(east, north), proj4string=CRS("+init=epsg:21097 +units=km"))
    lonLatCoords = spTransform(coordsUTM, CRS("+proj=longlat"))
    out = attr(lonLatCoords, "coords")
  }
  
  out
}

##### convert Binomial data to Bernoulli
# inputs:
# binDat: a binomial data frame
# outputs:
# each row of binDat is expanded into multiple Bernoulli trials.  Additional 
# variables are added to the data
binToBern = function(binDat, nVar="numChildren", yVar="died", weightVar="samplingWeight") {
  tryCatch({weights = binDat[,weightVar]})
  # dataForDirect.Rdata
}

##### simulate from a logistic-normal-binomial distribution
rLogisticNormBin = function(nsim, n=1, logitProb, logitProbVar=0) {
  # first make the parameter vectors the correct length (they should each be nsim long)
  makeRightLength = function(vals) {
    currentLength = length(vals)
    numReps = nsim/currentLength
    
    # make sure the user didn't mess up
    if(round(numReps) != numReps)
      stop("nsim not multiple of all lengths of argument vectors")
    else if(numReps < 1)
      stop("nsim has longer length than other arguments")
    
    rep(vals, numReps)
  }
  n = makeRightLength(n)
  logitProb = makeRightLength(logitProb)
  logitProbVar = makeRightLength(logitProbVar)
  
  # simulate matrix of all bernoulli values in each binomial random variable
  N = sum(n)
  bernVar = rep(logitProbVar, n)
  bernProb = expit(rep(logitProb, n) + rnorm(N, 0, sqrt(bernVar)))
  bernVals = rbinom(sum(n), rep(1, N), bernProb)
  
  # function for simulating individual logistic normal binomial random variable from bernoulli values
  binIDs = as.factor(rep(1:length(n), n))
  # aggregate(bernVals, list(binIDs=binIDs), sum)$x # apparently the tapply function is wayyyyy faster than aggregate...
  tapply(bernVals, binIDs, sum)
}

# parallel matrix multiply
# leftMat %*% rightMat
# cl: the cluster, as initialized via:
# cores=detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)
parMatMult = function(leftMat, rightMat, cl) {
  if (ncol(leftMat) != nrow(rightMat)) stop("Matrices do not conforme")
  idx   <- splitIndices(nrow(leftMat), length(cl))
  leftMatlist <- lapply(idx, function(ii) leftMat[ii,,drop=FALSE])
  ## ans   <- clusterApply(cl, leftMatlist, function(aa, rightMat) aa %*% rightMat, rightMat)
  ## Same as above, but faster:
  ans   <- clusterApply(cl, leftMatlist, get("%*%"), rightMat)
  do.call(rbind, ans)
}

# draw random numbers from an ecdf object
recdf = function(n, distribution) {
  probs = runif(n)
  quantile(distribution, probs, type=1) # type=1 signifies inverse ecdf
}

# get expectation of an ecdf object
ecdfExpectation = function(distribution) {
  distributionKnots = knots(distribution)
  distributionKnots = c(distributionKnots[1] - 1, distributionKnots)
  probs = distribution(distributionKnots[2:length(distributionKnots)]) - distribution(distributionKnots[1:(length(distributionKnots) - 1)])
  sum(distributionKnots[2:length(distributionKnots)] * probs)
}

# this function generates the sampling 
# probabilities of EAs within this strata
# i: the strata index, 1-47, corresponding to the row of easpc
# Nsamples: the total number of EAs drawn among strata with this urban/rural value
# m: the minimum number of draws in each strata
# urban: whether or not we are sampling from urban or rural strata
# log: return log probability or just probability
getStratumSamplingProb = function(i, nHHi, nHHTotal, Nsamples=617, m=3, urban=TRUE, dolog=FALSE) {
  # the number of households in this stratum (the commented out line should be the same as the one below)
  # Ni = ifelse(urban, easpc[i, "HHUrb"], easpc[i, "HHRur"])
  Ni = sum(nHHi)
  
  # the total number of households in Kenya that are either urban or rural
  Ntotal = ifelse(urban, sum(easpc$HHUrb), sum(easpc$HHRur))
  
  # the total number of strata with the given urban/rural categorization
  Nstrata = ifelse(urban, 47, 45)
  
  # the probability that you pick a particular ea in the first m draws from this stratum
  # pickEAinFirstM = dhyper(1, 1, Ni - 1, m)
  pickEAinFirstM = m * nHHi / Ni
  
  # the probability that you pick a particular ea in k additional draws from the stratum, 
  # given it wasn't drawn in the first m draws
  pickEAInK = function(k) {
    # successStates = Ni - m
    # totalStates = Ntotal - m * Nstrata
    # failureStates = totalStates - successStates
    totalDrawsLeft = Nsamples - m * Nstrata
    # probDrawStrataKTimes = dhyper(k, successStates, failureStates, totalDrawsLeft, log=TRUE)
    probDrawStrataKTimes = dbinom(k, totalDrawsLeft, Ni / Ntotal, log=TRUE)
    
    # countySuccessStates = 1
    # countyTotalStates = Ni - m
    # countyFailureStates = countyTotalStates - countySuccessStates
    # probDrawEaOutOfK = dhyper(1, countySuccessStates, countyFailureStates, k, log=TRUE)
    probDrawEaOutOfK = log(k * nHHi / Ni)
    
    # multiply the probabilities by summing on the log scale
    exp(probDrawStrataKTimes + probDrawEaOutOfK)
  }
  
  # the probability of picking this ea after the first m are drawn from the stratum 
  # is the probability it is drawn after the first m given that it was not drawn in 
  # the first m times the probability that it was not drawn in the first m
  maxDraws = min(Ni - m, Nsamples - m*Nstrata)
  pickEAafterFirstM = sum(sapply(1:maxDraws, pickEAInK)) * (1 - pickEAinFirstM)
  
  # probabilities of mutually exclusive events are summed
  pickEA = pickEAinFirstM + pickEAafterFirstM
  
  ifelse(dolog, log(pickEA), pickEA)
}

# this function generates the sampling probabilities of all households
# i: the strata index, 1-47, corresponding to the row of easpc
# Nsamples: the total number of EAs drawn among strata with this urban/rural value
# m: the minimum number of draws in each strata
# urban: whether or not we are sampling from urban or rural strata
# log: return log probability or just probability
getStratumSamplingProb2 = function(eaDat, eaDatLong, nSamplesUrban=617, nSamplesTotal=1612, m=3) {
  # the number of households in each stratum
  Nis = eaDat[, .(nHH=sum(nHH)), by=.(admin1, urban)]
  Ntotal = sum(Nis[["nHH"]])
  
  # the total number of households in Kenya that are either urban or rural
  Ntotal = ifelse(urban, sum(easpc$HHUrb), sum(easpc$HHRur))
  
  # the total number of strata with the given urban/rural categorization
  Nstrata = ifelse(urban, 47, 45)
  
  # the probability that you pick a particular ea in the first m draws from this stratum
  # pickEAinFirstM = dhyper(1, 1, Ni - 1, m)
  pickEAinFirstM = m * nHHi / Ni
  
  # the probability that you pick a particular ea in k additional draws from the stratum, 
  # given it wasn't drawn in the first m draws
  pickEAInK = function(k) {
    # successStates = Ni - m
    # totalStates = Ntotal - m * Nstrata
    # failureStates = totalStates - successStates
    totalDrawsLeft = Nsamples - m * Nstrata
    # probDrawStrataKTimes = dhyper(k, successStates, failureStates, totalDrawsLeft, log=TRUE)
    probDrawStrataKTimes = dbinom(k, totalDrawsLeft, Ni / Ntotal, log=TRUE)
    
    # countySuccessStates = 1
    # countyTotalStates = Ni - m
    # countyFailureStates = countyTotalStates - countySuccessStates
    # probDrawEaOutOfK = dhyper(1, countySuccessStates, countyFailureStates, k, log=TRUE)
    probDrawEaOutOfK = log(k * nHHi / Ni)
    
    # multiply the probabilities by summing on the log scale
    exp(probDrawStrataKTimes + probDrawEaOutOfK)
  }
  
  # the probability of picking this ea after the first m are drawn from the stratum 
  # is the probability it is drawn after the first m given that it was not drawn in 
  # the first m times the probability that it was not drawn in the first m
  maxDraws = min(Ni - m, Nsamples - m*Nstrata)
  pickEAafterFirstM = sum(sapply(1:maxDraws, pickEAInK)) * (1 - pickEAinFirstM)
  
  # probabilities of mutually exclusive events are summed
  pickEA = pickEAinFirstM + pickEAafterFirstM
  
  ifelse(dolog, log(pickEA), pickEA)
}

# get the probability of drawing a given household out of an enumeration area
getHHprob = function(nHH, nDraws=25, log=FALSE) {
  # dhyper(1, 1, nHH - 1, nDraws, log=log)
  nDraws/nHH
}

# match first occurrence of x in fun(values) using binary search assuming monotonicity of values
binarySearchMatch = function(x, values, fun=function(v, i) {v[i]}) {
  
  leftI = 1
  rightI = length(values)
  while(leftI < rightI) {
    matchI = floor((leftI + rightI) / 2)
    if(x <= fun(values, matchI))
      rightI = matchI
    else
      leftI = matchI + 1
  }
  
  if(x <= fun(values, leftI))
    leftI
  else
    NA
}

matchMultiple = function(vals, x) {
  x.map <- split(1:length(x), match(x, vals))
  # write a wrapper function that does a look-up on the unique list. and then returns all matches using the map.
  matchMake <- function(x) { x.map[[ fmatch(x, vals) ]] }
  sapply(vals, function(val) {matchMake(val)})
}

matchMultiple2 = function(vals, x) {
  dt_sol <- data.table(b, row = 1:length(b))[data.table(a), .(a, idx = row), on = .(b = a)][, .(idxs = list(idx)), by = a];
  dt_sol2 <- setNames(dt_sol$idxs, dt_sol$a)
  dt_sol2
}









