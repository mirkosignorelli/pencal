library(MASS) # for mvrnorm
library(Matrix) # for bdiag
source('functions/generate-t-weibull-v1.0.R')
source('functions/sim-eff-sizes.R')

# reference for name of dataset to save
simref = 4
name.part1 = paste('data/', simref, '/sim', simref, '_', sep = '')

# simulation settings
lambda = 0.1
nu = 2.2

n = 100
n.covs = 150
n.relev.covs = 10

set.seed(1234)
beta0 = runif(n.covs, min = 3, max = 7)
beta1 = sim.eff.sizes(n = n.covs, 
                          abs.range = c(1, 2))
eff.size.surv = sim.eff.sizes(n = n.relev.covs, 
                              abs.range = c(0.5, 1))
gamma = c(eff.size.surv, rep(0, n.covs - n.relev.covs))
delta = gamma

n.rep = 100
n.cases = 2 # few & many repeated measurements

for (a in 1:n.cases) {
  censoring.range = c(0.5, 10)
  if (a == 1) t.values = c(0, 0.5, 1, 2)
  if (a == 2) t.values = c(0, 0.1, 0.5, 1, 1.5, 2)
  for (h in 1:n.rep) {
    set.seed(h)
    # generate random effects
    Sigma = diag(1, n.covs)
    rand.int = mvrnorm(n = n, mu = rep(0, n.covs), Sigma = Sigma)
    rand.slope = mvrnorm(n = n, mu = rep(0, n.covs), Sigma = Sigma)
    m = length(t.values) # maximum theoretical number of repeated measurements
    id = rep(1:n, each = m)
    t.measures = rep(t.values, n)
    # generate protein values
    proteins = matrix(NA, n*m, n.covs)
    for (i in 1:n.covs){
      inter = rep(beta0[i] + rand.int[,i], each = m)
      slope = rep(beta1[i] + rand.slope[,i], each = m)
      mu = inter + slope*t.measures
      proteins[,i] = rnorm(n = n*m, mean = mu, sd = 0.7)
    }
    
    # simulate survival times
    X.temp = cbind(rand.int, rand.slope)
    true.t = generate.t.weibull(n = n, lambda = lambda, nu = nu,
            beta = c(gamma, delta), X = X.temp, seed = h)
    censoring.times = generate.censoring.times(n = n, 
                        min = censoring.range[1],
                        max = censoring.range[2])
    event = ifelse(true.t > censoring.times, 0, 1)
    observed.t = ifelse(event == 1, true.t, censoring.times)
    perc.censoring = length(which(event == 0)) / length(event)
  
    filename = paste(name.part1, a, '_', h, '.RData', sep = '')
    save(n, n.covs, n.relev.covs, proteins, 
         t.measures, t.values,
         m, id, rand.int, rand.slope,
         true.t, censoring.times, event, observed.t, perc.censoring,
         file = filename)
  }
  print(paste('case', a, 'done'))
}
