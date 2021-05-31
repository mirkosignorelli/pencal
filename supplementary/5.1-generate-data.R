library(MASS) # for mvrnorm
library(Matrix) # for bdiag
source('functions/generate-t-weibull-v1.0.R')
source('functions/sim-eff-sizes.R')

# reference for name of dataset to save
simref = 5
name.part1 = paste('data/', simref, '/sim', simref, '_', sep = '')

# simulation settings
lambda = 0.1
nu = 1.7

n = 500
# use first 300 subjects as training, the remaining 200 as validation
n.covs = 3
n.relev.covs = 2

set.seed(1234)
beta0 = c(4, 5, 6)
beta1 = c(1.5, -1, -1.2)
gamma = c(0.8, -0.6, 0)

n.rep = 100
n.cases = 1

for (a in 1:n.cases) {
  censoring.range = c(1.5, 10)
  if (a == 1) t.values = c(0, 0.1, 0.3, 0.5, 1)
  for (h in 1:n.rep) {
    set.seed(h)
    # generate random effects
    Sigma = diag(1, n.covs)
    rand.int = mvrnorm(n = n, mu = rep(0, n.covs), Sigma = Sigma)
    #rand.slope = mvrnorm(n = n, mu = rep(0, n.covs), Sigma = Sigma)
    m = length(t.values) # maximum theoretical number of repeated measurements
    id = rep(1:n, each = m)
    t.measures = rep(t.values, n)
    # generate protein values
    proteins = matrix(NA, n*m, n.covs)
    for (i in 1:n.covs){
      inter = rep(beta0[i] + rand.int[,i], each = m)
      #slope = rep(beta1[i] + rand.slope[,i], each = m)
      slope = rep(beta1[i], each = m)
      mu = inter + slope*t.measures
      proteins[,i] = rnorm(n = n*m, mean = mu, sd = 0.7)
    }
    
    # simulate survival times (from t = 1)
    true.t = 1 + generate.t.weibull(n = n, lambda = lambda, nu = nu,
            beta = gamma, X = rand.int, seed = h)
    #X.temp = cbind(rand.int, rand.slope)
    #true.t = generate.t.weibull(n = n, lambda = lambda, nu = nu,
    #        beta = c(gamma, delta), X = X.temp, seed = h)
    censoring.times = generate.censoring.times(n = n, 
                        min = censoring.range[1],
                        max = censoring.range[2])
    event = ifelse(true.t > censoring.times, 0, 1)
    observed.t = ifelse(event == 1, true.t, censoring.times)
    perc.censoring = length(which(event == 0)) / length(event)
    
    # visual checks:
    if (!on.shark & a == 1 & h == 1) {
      library(ptmixed)
      pdf(paste('figs-descr/', simref, '.1-traj.pdf'), 
          width = 10, height = 4)
      par(mfrow = c(1, 3))
      for (i in 1:3) {
        df = data.frame(x = t.measures, y = proteins[,i], id = id)
        make.spaghetti(x, y, id, data = df, ylim = c(-2, 10))
      }
      dev.off()
      library(corrplot)
      corrplot(cor(proteins))
      library(survival)
      kaplan <- survfit(Surv(time = observed.t, event = event) ~ 1,  
                        type="kaplan-meier")
      pdf(paste('figs-descr/', simref, '.1-KM.pdf'), 
          height = 6, width = 6)
      par(mar = c(4.5,5,2,2), bty = 'l')
      plot(kaplan, conf.int = F,
           xlab="t", ylab=expression(hat(S)(t)),
           col = 'blue', lwd = 1.3, cex.lab = 1.3, cex.axis = 1.3)
      abline(h = 0.2*c(0:5), lty = 2)
      dev.off()
    }
  
    filename = paste(name.part1, a, '_', h, '.RData', sep = '')
    save(n, n.covs, n.relev.covs, proteins, 
         t.measures, t.values,
         m, id, rand.int, #rand.slope,
         true.t, censoring.times, event, 
         observed.t, perc.censoring,
         file = filename)
  }
  print(paste('case', a, 'done'))
}

cens.distr = matrix(NA, n.rep, n.cases)
freq.1rmonly = mean.nrm = cens.distr

for (a in 1:n.cases) {
  for (h in 1:n.rep) {
    filename = paste(name.part1, a, '_', h, '.RData', sep = '')
    load(filename)
    # percentage of censored observations
    cens.distr[h, a] = perc.censoring
    min.c.t = apply(cbind(true.t, censoring.times), 1, min)
    # how many individuals have only 1 RM
    freq.1rmonly[h, a] = sum(min.c.t <= t.values[2])
    # mean number of RMs
    n.rm = rep(NA, n)
    for (i in 1:n) n.rm[i] = sum(t.values < min.c.t[i])
    mean.nrm[h, a] = mean(n.rm)
  }
}

# censoring distribution
( avg.cens = round(100*apply(cens.distr, 2, mean),0) )
# mean number of RMs
( avg.nrm = round(apply(mean.nrm, 2, mean),1) )
# % of subjects with just 1 RM
avg.1rmonly = round(apply(freq.1rmonly, 2, mean),0)/n*100
round(avg.1rmonly, 1)

save(avg.cens, avg.nrm, avg.1rmonly,
  file = paste('data/features_sim', simref, '.RData', sep = ''))
