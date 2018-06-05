library(haven)
library(survey)

# -- a few helpful functions -- #
expit<-function(x){
  exp(x)/(1+exp(x))
}

logit<-function(x){
  log(x/(1-x))
}

# function to extend dataset to binary form 
extendData <- function(clustDatRow, v001){
  
  # add extra columns for ageMonth, ageGrpD, v001, v002
  nC = clustDatRow$numChildren
  tmp = data.frame(clustDatRow[c(1, 6:16)])
  tmp$v001 = v001

  ageMonth = rep(0, nC)
  ageGrpD = rep("[0,1)", nC)
  v001 = rep(v001, nC)
  # there is only one child and one mother per household.
  # All 25 households are sampled
  v002 = 1:nC
  
  died = c(rep(0,nC-clustDatRow$died), rep(1, clustDatRow$died))
  if(clustDatRow["urban"][1,1]){
    urbanRural = rep("urban", nC)
  } else {
    urbanRural = rep("rural", nC)
  }

  res = merge(data.frame(died, ageMonth, ageGrpD, v001, v002, urbanRural), tmp, by="v001")
  
  res$regionRural <- with(res, interaction(admin1, urbanRural), drop=TRUE)
  
  return(res)
}


# - a function that reads in a glm or svyglm - #
# - object and returns the estimate and SE - #
# - specifics in the supplementary materials - #
## This function takes care of the delta method
## to calculate the variance of u5m as a function
## of the age specific hazards, \beta_a .

get.est<-function(glm.ob){

  beta<-summary(glm.ob)$coef[,1]
 
  u1m.est <-expit(beta)
  var.est <- vcov(glm.ob)[1,1]
  
  # compute 80% CI intervals
  lower <- logit(u1m.est)+qnorm(c(0.9))*sqrt(var.est)
  upper <- logit(u1m.est)+qnorm(c(0.1))*sqrt(var.est)
  return(c(u1m.est,lower, upper,logit(u1m.est),var.est))
}

# -- a function to subset the design based on a region and time period -- #
# -- and then run the svyglm function and return the get.est() results -- #

## First line in function allows you to subset your data and ALSO the specified
## svydesign object into area (usually v024 variable in DHS) 
## and time (per5 is a variable we construct for the 5-year periods in the Stata step)
## Second line fits the survey-weighted glm

region.time.HT<-function(dataobj, svydesign, area){
  
  tmp<-subset(svydesign, (admin1==area))
  
  tt2 <- tryCatch(glmob<-svyglm(died.x~1,
                                design=tmp,family=quasibinomial, maxit=50), 
                  error=function(e) e, warning=function(w) w)
  
  if(is(tt2, "warning")){
    if(grepl("agegroups", tt2)){
      res <- get.est(glmob)
      res = c(res, 2)
    } else {
      res = c(rep(NA, 5), 3)
    }
    return(res)
  }
  if(is(tt2,"error")){
    res = c(rep(NA, 5), 1)
    return(res)
  } else {
    res <- get.est(glmob)
    res = c(res, 0)
    return(res)
  }
}



defineSurvey <- function(childBirths_obj, stratVar, useSamplingWeights=TRUE){
  
  options(survey.lonely.psu="adjust")
  
  # --- setting up a place to store results --- #
  regions <- sort(unique(childBirths_obj$admin1))
  regions_num  <- 1:length(regions)
  
  results<-data.frame(admin1=rep(regions,each=1))
  results$var.est<-results$logit.est<-results$upper<-results$lower<-results$u1m<-NA
  results$converge <- NA
  
  if(useSamplingWeights){
    childBirths_obj$wt <- childBirths_obj$samplingWeight
  } else {
    childBirths_obj$wt <- NULL
  }

  if(is.null(stratVar)){
    # --- setting up the design object --- #
    ## NOTE: -the v001 denote
    ##        one stage cluster design (v001 is cluster)
    ##       -This call below specifies our survey design
    ##        nest = T argument nests clusters within strata
    my.svydesign <- svydesign(id= ~v001,
                              strata =NULL,
                              weights=NULL, data=childBirths_obj)
  } else {
    ## not in all surveys does v022 contain the correct sampling strata
    ## Thus, the correct vector has to be provided externally
    childBirths_obj$strat <- stratVar
  
    # --- setting up the design object --- #
    ## NOTE: -the v001 denote
    ##        one stage cluster design (v001 is cluster)
    ##       -This call below specifies our survey design
    ##        nest = T argument nests clusters within strata
    my.svydesign <- svydesign(id= ~v001,
                              strata=~strat, nest=T, 
                              weights=~wt, data=childBirths_obj)
  }
  
  for(i in 1:nrow(results)){
    results[i, 2:7] <- region.time.HT(dataobj=childBirths_obj, svydesign=my.svydesign, 
                                      area=results$admin1[i])
  }
  return(results)
}


run_naive <- function(childBirths_obj){
  
  regions <- sort(unique(childBirths_obj$admin1))
  regions_num  <- 1:length(regions)
  
  results<-data.frame(admin1=rep(regions,each=1))
  results$var.est<-results$logit.est<-results$upper<-results$lower<-results$u1m<-NA
  results$converge <- NA
  
  for(i in 1:nrow(results)){
    my.glm <- glm(died.x~1, family=binomial, 
                  data=childBirths_obj, 
                  subset = admin1 == results$admin1[i] ) 
    # newdat = childBirths_obj[childBirths_obj$admin1==results$admin1[i], ]
    # my.glm2 <- glm(died.x~1, family=binomial, 
    #               data=newdat) 
    
    results[i, 2:7] <- c(get.est(my.glm),0)
  }
  return(results)
}

