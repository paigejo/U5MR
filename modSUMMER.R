myCountrySummary = function(births, years, idVar = "v002", regionVar = "region", 
          timeVar = "per5", clusterVar = "~v001+v002", ageVar = "ageGrpD", 
          weightsVar = "v005", geo.recode = NULL) 
{
  if (is.null(births)) {
    stop("No births file specified!")
  }
  if (is.null(years)) {
    stop("Names of the 5-year intervals not provided!")
  }
  if (!regionVar %in% colnames(births)) {
    if ("v101" %in% colnames(births)) {
      colnames(births)[which(colnames(births) == "v101")] <- regionVar
      warning("region variable not defined: using v101", 
              immediate. = TRUE)
    }
    else if ("v024" %in% colnames(births)) {
      colnames(births)[which(colnames(births) == "v024")] <- regionVar
      warning("region variable not defined: using v024", 
              immediate. = TRUE)
    }
    else {
      stop("region variable not defined, and no v101 or v024!")
    }
  }
  if (!is.null(geo.recode)) {
    births <- ChangeRegion(births, Bmat = geo.recode, regionVar = regionVar)
  }
  births$region0 <- births[, regionVar]
  births$id0 <- births[, idVar]
  births$weights0 <- births[, weightsVar]
  births$time0 <- births[, timeVar]
  births$age0 <- births[, ageVar]
  if (sum(c(0, 1, 12, 24, 36, 48) %in% births$age0) == 6) {
    births$age0[births$age0 == 0] <- "0"
    births$age0[births$age0 == 1] <- "1-11"
    births$age0[births$age0 == 12] <- "12-23"
    births$age0[births$age0 == 24] <- "24-35"
    births$age0[births$age0 == 36] <- "36-47"
    births$age0[births$age0 == 48] <- "48-59"
  }
  time_inconsistent <- which(!(births$time0 %in% years))
  if (length(time_inconsistent) > 0) {
    warning(paste("Name for time periods are inconsistent. Found the following levels in data:", 
                  unique(births$time0[time_inconsistent])), immediate. = TRUE)
  }
  if (is.null(births$strata)) {
    stop("Strata not defined.")
  }
  options(survey.lonely.psu = "adjust")
  my.svydesign <- survey::svydesign(ids = ~id0, cluster = stats::formula(clusterVar), 
                                    strata = ~strata, nest = T, weights = ~weights0, data = births)
  regions_list <- as.character(sort(names(table(births$region0))[as.vector(table(births$region0) != 
                                                                             0)]))
  regions_num <- 1:length(regions_list)
  regions_list <- c("All", regions_list)
  regions_num <- c(0, regions_num)
  results <- data.frame(region = rep(regions_list, each = length(years)))
  results$region_num <- rep(regions_num, each = length(years))
  results$years <- rep(years, length(regions_list))
  results$var.est <- results$logit.est <- results$upper <- results$lower <- results$u5m <- NA
  region.time.HT.withNA <- function(which.area, which.time) {
    time0 <- NULL
    region0 <- NULL
    if (which.area == "All") {
      tmp <- subset(my.svydesign, (time0 == which.time))
    }
    else {
      tmp <- subset(my.svydesign, (time0 == which.time & 
                                     region0 == as.character(which.area)))
    }
    if (dim(tmp)[1] == 0) {
      return(rep(NA, 5))
    }
    else if (sum(tmp$variables$died) == 0) {
      warning(paste0(which.area, " ", which.time, " has no death, set to NA\n"), 
              immediate. = TRUE)
      return(rep(NA, 5))
    }
    else {
      glm.ob <- survey::svyglm(died ~ 1, 
                               design = tmp, family = stats::quasibinomial, 
                               maxit = 50)
      return(get.est.withNA(glm.ob))
    }
  }
  get.est.withNA <- function(glm.ob) {
    V <- matrix(0, 1, 1)
    betas <- rep(NA, 1)
    # labels <- c("0", "1-11", "12-23", "24-35", "36-47", "48-59")
    labels = c("(Intercept)")
    # labels <- paste("factor(age0)", labels, sep = "")
    colnames(V) <- rownames(V) <- labels
    names(betas) <- labels
    V2 <- stats::vcov(glm.ob)
    if (length(which(colnames(V2) %in% colnames(V) == FALSE)) > 
        0) {
      stop("Error for input age group names!")
    }
    V[rownames(V2), colnames(V2)] <- V2
    betas2 <- summary(glm.ob)$coef[, 1]
    # betas[names(betas2)] <- betas2
    betas = betas2
    # ns <- c(1, 11, 12, 12, 12, 12)
    ns <- c(1)
    probs <- expit(betas)
    u5m.est <- (1 - prod((1 - probs)^ns, na.rm = TRUE))
    gamma <- prod((1 + exp(betas))^ns, na.rm = TRUE)
    derivatives <- (gamma)/(gamma - 1) * ns * expit(betas)
    derivatives[which(is.na(derivatives))] <- 0
    var.est <- t(derivatives) %*% V %*% derivatives
    lims <- logit(u5m.est) + stats::qnorm(c(0.025, 0.975)) * 
      sqrt(c(var.est))
    return(c(u5m.est, expit(lims), logit(u5m.est), var.est))
  }
  x <- mapply(region.time.HT.withNA, which.area = results$region, 
              which.time = results$years)
  results[, 4:8] <- t(x)
  return(results)
}