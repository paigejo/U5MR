# load direct estimates
# load("resultsDirectNaive.RData")
# load a different 1 of these depending on whether a cluster effect should be included 
# in the simulation of the data or not (tausq is the cluster effect variance)
# load("resultsDirectNaiveTausq0.RData")
load("resultsDirectNaiveTausq0.01.RData")

# source Mercer et al. code
source("mercer.R")

# number of simulation scenarios
n = 100

# go throught all of them 
library(INLA)

mercerSRS = merceroverSamp = list()

for(i in 1:100){
  print(i)
  tmpoverSamp = mercer_u1m(directEstoverSamp[[i]]$logit.est, directEstoverSamp[[i]]$var.est, 
                graph.path = "Kenyaadm1.graph")
  resoverSamp=   data.frame(admin1=directEstoverSamp[[i]]$admin1, 
                       u1m.mercer=expit(tmpoverSamp$summary.linear.predictor$mean),
                       lower.mercer=tmpoverSamp$summary.linear.predictor$"0.1quant",
                       upper.mercer=tmpoverSamp$summary.linear.predictor$"0.9quant",
                       logit.est.mercer=tmpoverSamp$summary.linear.predictor$mean, 
                      var.est.mercer=(tmpoverSamp$summary.linear.predictor$sd)^2)
  merceroverSamp[[i]] = resoverSamp
  
  
  tmpSRS = mercer_u1m(directEstSRS[[i]]$logit.est, directEstSRS[[i]]$var.est, 
                    graph.path = "Kenyaadm1.graph")
  
  resSRS = data.frame(admin1=directEstSRS[[i]]$admin1,
                  u1m.mercer=expit(tmpSRS$summary.linear.predictor$mean),
                  lower.mercer=tmpSRS$summary.linear.predictor$"0.1quant",
                  upper.mercer=tmpSRS$summary.linear.predictor$"0.9quant",
                  logit.est.mercer=tmpSRS$summary.linear.predictor$mean, 
                  var.est.mercer=(tmpSRS$summary.linear.predictor$sd)^2)
  mercerSRS[[i]] = resSRS
}

# save(merceroverSamp, mercerSRS, file="resultsMercer.RData")
# save(merceroverSamp, mercerSRS, file="resultsMercerTausq0.RData")
save(merceroverSamp, mercerSRS, file="resultsMercerTausq0.01.RData")
