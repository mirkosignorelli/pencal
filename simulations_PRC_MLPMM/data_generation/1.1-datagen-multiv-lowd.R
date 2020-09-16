library(MASS) # for mvrnorm
library(Matrix) # for bdiag
source('functions/generate-t-weibull-v1.0.R')
source('functions/sim-eff-sizes.R')

# reference for name of dataset to save
simref = 1
name.part1 = paste('data/', simref, '/sim', simref, '_', sep = '')

# simulation settings
lambda = 0.1
nu = 1.7

n = 100
n.prot = 10
n.antib = 3
n.relev.prot = 4

set.seed(1234)
beta0 = runif(n.prot*n.antib, min = 3, max = 7)
beta1 = sim.eff.sizes(n = n.prot*n.antib, 
                          abs.range = c(1, 2))
eff.size.surv = sim.eff.sizes(n = n.relev.prot, 
                              abs.range = c(0.5, 1))
gamma = c(eff.size.surv, rep(0, n.prot - n.relev.prot))
delta = gamma

n.rep = 100
n.cases = 1

for (a in 1:n.cases) {
  censoring.range = c(0.5, 10)
  if (a == 1) t.values = c(0, 0.1, 0.5, 1, 1.5, 2)
  for (h in 1:n.rep) {
    set.seed(h)
    # generate random effects
    sigma.b0 = diag(0.5, n.prot*n.antib)
    b0 = mvrnorm(n = n, mu = rep(0, n.prot*n.antib), Sigma = sigma.b0)
    sigma.u = matrix(c(1, 0.5, 0.5, 1), 2, 2)
    sigma = sigma.u
    for (i in 2:n.prot) sigma = bdiag(sigma, sigma.u)
    us = mvrnorm(n = n, mu = rep(0, n.prot*2), Sigma = sigma)
    u0 = us[, seq(1, 2*n.prot-1, by = 2)]
    u1 = us[, seq(2, 2*n.prot, by = 2)]
    m = length(t.values) # maximum theoretical number of repeated measurements
    id = rep(1:n, each = m)
    t.measures = rep(t.values, n)
    # generate protein values
    markers = matrix(NA, n*m, n.prot*n.antib)
    for (i in 1:(n.prot*n.antib)){
      j = ceiling(i / 3)
      inter = rep(beta0[i] + u0[,j] + b0[,i], each = m)
      slope = rep(beta1[i] + + u1[,j], each = m)
      mu = inter + slope*t.measures
      markers[,i] = rnorm(n = n*m, mean = mu, sd = 0.7)
    }
    
    # simulate survival times
    X.temp = cbind(u0, u1)
    true.t = generate.t.weibull(n = n, lambda = lambda, nu = nu,
            beta = c(gamma, delta), X = X.temp, seed = h)
    censoring.times = generate.censoring.times(n = n, 
                        min = censoring.range[1],
                        max = censoring.range[2])
    event = ifelse(true.t > censoring.times, 0, 1)
    observed.t = ifelse(event == 1, true.t, censoring.times)
    perc.censoring = length(which(event == 0)) / length(event)

    filename = paste(name.part1, a, '_', h, '.RData', sep = '')
    save(n, n.prot, n.antib, n.relev.prot, markers, 
         t.measures, t.values,
         m, id, b0, u0, u1,
         true.t, censoring.times, event, observed.t, perc.censoring,
         file = filename)
  }
  print(paste('case', a, 'done'))
}