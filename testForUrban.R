

# extract births in the range 2005 to 2014 (most recent five years)
cutOff = 0
cutOff = 11
lows = seq(1980, 2010, by=5)
lows = seq(1979, 2009, by=5)
lows = seq(1980, 2009, by=5)
lows = seq(1978, 2008, by=5)
for(i in 1:length(lows)) {
  subdata <- data.frame(data[,c('b2', 'b1', 'b5', 'b7', 'v024', 'v025', 'v001', 'v002', 'v005', 'v006', 'v007')])
  lowYear<-lows[i]
  highYear<-lowYear + 4
  subdata <- subdata[(subdata[,'b2'] >= lowYear & subdata[,'b2'] <= highYear),]
  
  # add a column for the stratification variable as an interaction between
  # the urban/rural indicator 'v025' (1: urban, 2:rural) and the region indicator 'v024'
  subdata$regionUral <- with(subdata, interaction(v024, v025), drop=TRUE)
  
  # add a column for the unique households with interaction between
  # the household indicator 'v002' and the cluster indicator 'v001'
  subdata$hhold <- with(subdata, interaction(v001, v002), drop=TRUE)
  
  # find for each cluster the regionUral indicator
  clStrat = subdata[,c("v001", "regionUral", "v024", "v025", "v005")]
  clStrat = clStrat[!duplicated(clStrat), ]
  colnames(clStrat) = c("clusterID", "regionRural", "region", "urban", "samplingWeight")
  clStrat$urban =  clStrat$urban == 1
  
  # determine whether each child survived for at least one month
  lived = subdata$b7 > cutOff
  lived[is.na(lived)] = TRUE
  died = !lived
  
  # get the number of birth by cluster
  n <- table(subdata[,'v001'])
  clusterid <- dimnames(n)[[1]]
  n.data = data.frame(clusterID=clusterid, n=as.vector(n))
  
  # get the number of deaths by cluster
  # y <- table(subdata[,'v001'])
  y = aggregate(died, list(subdata[,'v001']), sum)
  n.data$y = y$x
  
  # add in strata
  mort <- merge(n.data, clStrat, by='clusterID', all=TRUE, sort=TRUE)
  
  # Read geographical information
  library(rgdal)
  spObj = readOGR(dsn = "Kenya2014gps/", layer = "KEGE71FL")
  
  # Extract (lon, lat) coordinates of all clusters
  geoObj = data.frame(cId = spObj$DHSCLUST, lon = spObj$LONGNUM, lat = spObj$LATNUM)
  
  # Extract coordinates of clusters with data
  idx = match(mort$clusterID, geoObj$cId)
  mort$lon = geoObj$lon[idx]
  mort$lat = geoObj$lat[idx]
  
  # Missing geographical information is assigned value (0,0)
  # Remove these
  missIdx = which(mort$lon == 0)
  mort = mort[-missIdx,]
  
  library(SUMMER)
  library(foreign)
  
  gpsDat = readShapePoints("Kenya2014gps/KEGE71FL.shp")
  coords = attr(gpsDat, "coords")
  plot(coords)
  world(add=TRUE)
  test = coords[,1] < 20 #  these are the observations  whose source is missing.  remove these
  sum(test)
  
  #  remove observations with unknown locations and set unknown data to NA
  gpsDat=gpsDat[!test,]
  names(gpsDat)=c("ID","countryID", "year", "clustID", "countryIDFIPS", "countryIDAdminID", 
                  "admin1FIPS","admin1IDSALB", 
                  "admin1SALB", "admin1ID", "admin1", "regionID", "region", 
                  "source", "urban", "lat", "lon", "altGPS", "altRadar", "coordRef")
  gpsDat$altGPS[gpsDat$altGPS == 9999] = NA
  gpsDat$altGPS[gpsDat$altRadar == 9999] = NA
  gpsDat$urban = gpsDat$urban == "U"
  gpsDat$countryIDFIPS[gpsDat$countryIDFIPS == "NULL"] = NA
  gpsDat$admin1FIPS[gpsDat$admin1FIPS == "NULL"] = NA
  gpsDat$admin1IDSALB[gpsDat$admin1IDSALB == "NULL"] = NA
  gpsDat$admin1SALB[gpsDat$admin1SALB == "NULL"] = NA
  gpsDat$countryIDAdminID[gpsDat$countryIDAdminID == "NULL"] = NA
  
  # get region and admin data from gps data, add to clusters in mort dataset
  gpsI = match(data.frame(rbind(mort$lon, mort$lat)), data.frame(rbind(gpsDat$lon, gpsDat$lat)))
  mort$admin1 = gpsDat$admin1[gpsI]
  mort$region = gpsDat$region[gpsI]
  
  # get easting and northing using projection
  tmp = projKenya(mort$lon, mort$lat)
  mort$east = tmp[,1]
  mort$north = tmp[,2]
  
  out = runBYM2Mort(mort, includeUrbanRural = TRUE, includeCluster = TRUE, saveResults = FALSE)
  if(i == 1) {
    parameters = list(out$parameters)
  } else {
    parameters = c(parameters, list(out$parameters))
  }
}
names(parameters) = lows
for(i in 1:length(lows))
  names(parameters)[i] = paste0(names(parameters)[i], "-", lows[i]+4)
test1m = lapply(parameters, function(x) {round(as.matrix(x), 3)})
test1y = lapply(parameters, function(x) {round(as.matrix(x), 3)})
test1y2 = lapply(parameters, function(x) {round(as.matrix(x), 3)})


