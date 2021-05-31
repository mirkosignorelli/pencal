h = 1 # 1 to 100 (parallelized on a HPC cluster)

library(pencal)
library(survival)
library(joineRML)
library(survivalROC)
library(foreach)

simref = 5
n.cases = 1

a = 1 
name.part1 = paste('data/', simref, '/sim', simref, '_', sep = '')
filename = paste(name.part1, a, '_', h, '.RData', sep = '')
load(filename)

### validation set: ###
n.training = 300
# datasets to supply to pencal:
id.valid = which(id > n.training)
long.data.valid = data.frame(id = id[id.valid], 
                       time = t.measures[id.valid], 
                       proteins[id.valid, ])
names(long.data.valid)[3:5] = c('X1', 'X2', 'X3')
surv.data.valid = data.frame(id = (n.training+1):n, 
                       time = observed.t[-c(1:n.training)],
                       event = event[-c(1:n.training)])

# dataset to supply to joineRML:
nreps.id = table(long.data.valid$id)
all.data.valid = data.frame(long.data.valid, 
             observed.t = rep(observed.t[-c(1:n.training)], times = nreps.id),
             event = rep(event[-c(1:n.training)], times = nreps.id))

subf = paste(simref, '.2', sep = '')
name.part2 = paste('results/',subf, '/res', subf, '_', sep = '')
resfile = paste(name.part2, a, '_', h, '.RData', sep = '')
load(resfile)

t.pred = 2:5

### predict conditional survival probabilities ###
# pencal with penalty
pred.uncond = survpred_prclmm(step1, step2, step3.ridge, 
                      times = t.pred,
                      new.longdata = long.data.valid,
                      new.basecovs = NULL, keep.ranef = T)
pred.t1 = survpred_prclmm(step1, step2, step3.ridge,
                      times = 1,
                      new.longdata = long.data.valid,
                      new.basecovs = NULL, keep.ranef = T)
pred.pencal.ridge = pred.uncond$predicted_survival
for (i in 2:ncol(pred.pencal.ridge)) {
  pred.pencal.ridge[ , i] = pred.uncond$predicted_survival[ , i] / pred.t1$predicted_survival[ , 2]
}
names(pred.pencal.ridge) = c('id', paste('S(', t.pred, '|1)', sep = ''))

# PRC without penalization:
df.v = data.frame(surv.data.valid, pred.uncond$predicted_ranefs)
temp.sfit = survfit(step3.coxph, newdata = df.v, 
                    se.fit = F, conf.int = F)
shat.uncond = t(summary(temp.sfit, times = t.pred)$surv)
shat.t1 = t(summary(temp.sfit, times = 1)$surv)
pred.pencal.nopenalty = shat.uncond

for (i in 1:ncol(shat.uncond)) {
  pred.pencal.nopenalty[ , i] = shat.uncond[ , i] / shat.t1
}
pred.pencal.nopenalty = data.frame(id = surv.data.valid$id, 
                                   pred.pencal.nopenalty)
names(pred.pencal.nopenalty) = c('id', paste('S(', t.pred, '|1)', sep = ''))

# joineRML
pred.joineRML = foreach(i = (n.training+1):n, .combine = 'rbind') %do% {
  df.sbj = subset(all.data.valid, id == i)
  p.hat = dynSurv(fit.joineRML, newdata = df.sbj,
                  u = t.pred)
  p.hat$pred$surv
}

pred.joineRML = data.frame((n.training+1):n, pred.joineRML)
names(pred.joineRML) = c('id', paste('S(', t.pred, '|1)', sep = ''))

### performance measures ###
t.pred = 2:5
tdauc.vals = data.frame(t.pred = t.pred,
                        joineRML = NA,
                        pencal.ridge = NA,
                        pencal.nopenalty = NA)

### compute time-dependent AUC ###
for (t in t.pred) {
  tdauc.vals[t-1, 2] = survivalROC(Stime = surv.data.valid$time, 
                               status = surv.data.valid$event, 
                               marker = 1 - pred.joineRML[ , t], 
                               predict.time = t, 
                               method = "NNE", span = 0.25*n^(-0.2))$AUC
  tdauc.vals[t-1, 3] = survivalROC(Stime = surv.data.valid$time, 
                             status = surv.data.valid$event, 
                             marker = 1 - pred.pencal.ridge[ , t], 
                             predict.time = t, 
                             method = "NNE", span = 0.25*n^(-0.2))$AUC
  tdauc.vals[t-1, 4] = survivalROC(Stime = surv.data.valid$time, 
                             status = surv.data.valid$event, 
                             marker = 1 - pred.pencal.nopenalty[ , t], 
                             predict.time = t, 
                             method = "NNE", span = 0.25*n^(-0.2))$AUC
}

subf = paste(simref, '.4', sep = '')
name.part3 = paste('results/',subf, '/res', subf, '_', sep = '')
resfile = paste(name.part3, a, '_', h, '.RData', sep = '')
save(pred.pencal.ridge, pred.joineRML,
     tdauc.vals, file = resfile)
