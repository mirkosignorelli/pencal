h = 1 # 1 to 100 (parallelized on a HPC cluster)

library(pencal)
library(survival)
library(joineRML)

simref = 5
n.cases = 1
name.part1 = paste('data/', simref, '/sim', simref, '_', sep = '')
subf = paste(simref, '.2', sep = '')
name.part2 = paste('results/',subf, '/res', subf, '_', sep = '')
a = 1 # no need to repeat a = 1, 2 for this model (same results...)

filename = paste(name.part1, a, '_', h, '.RData', sep = '')
load(filename)

# datasets to supply to pencal:
n.training = 300
id.train = which(id <= n.training)
long.data = data.frame(id = id[id.train], 
                       time = t.measures[id.train], 
                       proteins[id.train, ])
names(long.data)[3:5] = c('X1', 'X2', 'X3')
surv.data = data.frame(id = 1:n.training, 
                       time = observed.t[1:n.training],
                       event = event[1:n.training])

# pencal
t0 = Sys.time()
step1 = fit_lmms(y.names = c('X1', 'X2', 'X3'),
                 fixefs = ~ time,
                 ranefs = ~ 1 | id,
                 long.data = long.data, 
                 surv.data = surv.data,
                 t.from.base = time)
step2 = summarize_lmms(step1)
t1 = Sys.time()
# with ridge penalty in step 3:
step3.ridge = fit_prclmm(step2, surv.data = surv.data)
t2a = Sys.time()
# without penalization in step 3:
df.train.coxph = data.frame(surv.data, step2$ranef.orig)
x.names = names(step2$ranef.orig)
formula = as.formula(paste('Surv(time, event) ~', 
                paste(x.names, collapse = '+')))
t2b = Sys.time()
step3.coxph = coxph(formula, data = df.train.coxph)
t2c = Sys.time()

# dataset to supply to joineRML:
nreps.id = table(step1$df.sanitized$id)
all.data = data.frame(step1$df.sanitized[ , -6], 
      observed.t = rep(observed.t[1:n.training], times = nreps.id),
      event = rep(event[1:n.training], times = nreps.id))

t3 = Sys.time()
fit.joineRML = mjoint(
  formLongFixed = list("X1" = X1 ~ time, "X2" = X2 ~ time,
                       "X3" = X3 ~ time),
  formLongRandom = list("X1" = ~ 1 | id, "X2" = ~ 1 | id,
                        "X3" = ~ 1 | id),
  formSurv = Surv(observed.t, event) ~ 1,
  data = list(all.data, all.data, all.data),
  timeVar = "time")
t4 = Sys.time()

# record computing time:
timediff = function(t1, t2) {
  as.numeric(difftime(t2, t1, units = 'secs'))
}
time.eval = data.frame(method = NA, time = NA)
time.eval[1,] = c('pencal.ridge',
                  timediff(t0, t1)+timediff(t1, t2a))
time.eval[2,] = c('pencal.no.penalty',
                  timediff(t0, t1)+timediff(t2b, t2c))
time.eval[3,] = c('joineRML',
                  timediff(t3, t4))

resfile = paste(name.part2, a, '_', h, '.RData', sep = '')
save(step1, step2, step3.ridge, step3.coxph, fit.joineRML,
     long.data, surv.data, df.train.coxph, all.data, time.eval,
     file = resfile)

