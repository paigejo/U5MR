# get PC prior for the BYM2 total variance
library(INLA)
set.seed(123)
N = 10000
U = 1
prec = inla.pc.rprec(alpha = 0.01, u = U, n = N)
sig = 1/sqrt(prec)
z = rnorm(N, mean = rep(0, N), sd = sig)
effect = exp(z)
quantile(effect, c(0.025, 0.975))
#      2.5%     97.5% 
# 0.5016154 1.8766063 

# get PC prior for the cluster effect variance
Mortality = 0.04 # Kenya's first year mortality rate is approximately 4% for 2005-2009
mu = log(Mortality/(1-Mortality))
N = 1000000
U = 3
alpha = 0.01
prec = inla.pc.rprec(alpha = alpha, u = U, n = N)
sig = 1/sqrt(prec)
z = rnorm(N, mean = rep(mu, N), sd = sig)
p = 1/(1+exp(-z))
y = matrix(rbinom(n = 2*N, size = 1, prob = rep(p,each = 2)), ncol = 2, byrow = T)
print(cor(y)[1,2])
# 0.1130728