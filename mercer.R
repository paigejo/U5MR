# CAUTION: logit.est and var.est must be provided for admin1 in alphabetic order
# otherwise the graph file will not fit.
mercer_u1m = function(logit.est, var.est, graph.path){
  
  # our weighted estimates represent the outcome 
  # assumed to be normally distributed with the true mean and
  # fixed variance equal to var.est
  y = logit.est
  
  # our model includes an intercept, iid space and ICAR component
  # use hyperpriors as in Mercer et al. 
  formula = y ~ f(idx, model="iid", 
                  hyper=list(prec=list(param=c(0.5, 0.001488), prior="loggamma"))) + 
                f(idx2, model="besag", graph="Kenyaadm1.graph", 
                  hyper=list(prec=list(param=c(0.5, 0.00360), prior="loggamma")))
  
  ## generate dataset
  data = list(y=y, scale=1/var.est, idx=1:length(y), idx2=1:length(y))
  
  # be careful to fix the likelihood variance, since we have observation-specific
  # variances we need to use the argument scale and fix the overall precision to 1.
  result = inla(formula, family="gaussian", data=data, 
                control.family=list(hyper=list(prec=list(initial=0, fixed=TRUE))), 
                control.predictor=list(compute=TRUE, link=1, quantiles=c(0.025, 0.1, 0.5, 0.9, 0.975)),
                scale = scale)
}

