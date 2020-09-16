###########################
nsub <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(R.version.string)
print(paste('node(s):', Sys.getenv("SLURM_JOB_NODELIST")))
cat('simulation number', nsub, '\n') # simulation number
h = nsub
###########################

on.shark = T

if (on.shark) setwd('/exports/molepi/ms/hdorc2-ridge')
if (on.shark == F) {
  rm(list = ls())
  library(msigno)
  clear.screen(); set.cf()
  h = 1
  a = 2
  boot = 1
}

library(glmnet)
library(survival)
library(survcomp) # C index (concordance.index)
library(survivalROC) # time-dependent AUC
# load function to compute predicted random effects
# on test set using the model estimated on the training set:
source('./functions/ranef-predictors-v2.4.R')

# standardize predictors in glmnet?
do.stand = T

simref = 1
n.cases = 2
name.part1 = paste('data/', simref, '/sim', simref, '_', sep = '')
subf = paste(simref, '.4', sep = '')
name.part2 = paste('results/',subf, '/res', subf, '_', sep = '')

# number of bootstrap samples
n.boots = 100
# time(s) at which to evaluate tdAUC
t.pred = c(3, 5, 7)

for (a in 1:n.cases) {
  start.time = Sys.time()
  
  filename = paste(name.part1, a, '_', h, '.RData', sep = '')
  load(filename)
  survdata = Surv(time = observed.t, event = event)
  # remove the measurements a) taken after loss of ambulation
  # and b) after censoring
  temp = cbind(censoring.times, true.t)
  t.last = apply(temp, 1, min)
  t.last.rep = rep(t.last, each = m)
  rows.rem = which(t.measures > t.last.rep)
  id = id[- rows.rem]
  t.measures = t.measures[- rows.rem]
  proteins = proteins[-rows.rem, ]
  
  out = as.data.frame(matrix(NA, n.boots+1, 2+2+2*length(t.pred)))
  names(out) = c('repetition', 'boot', 'c.boot', 'c.orig',
                 paste(rep(c('auc.boot.t', 'auc.orig.t'), length(t.pred)), 
                       rep(t.pred, each = 2), sep = ''))
  
  out$repetition = h
  k = 1
  # begin bootstrap!  
  for (boot in 0:n.boots) {
    out$boot[k] = boot
    set.seed(boot)
    # boot = 0: estimate model on original dataset
    if (boot == 0) boot.ids = 1:n
    # boot > 0: actual bootstrap
    # sample individuals for bootstrap:
    if (boot > 0) {
      boot.ids = sample(1:n, n, replace = T)
      boot.ids = boot.ids[order(boot.ids)]
    }
    # prepare "boot" datasets on which the model will be trained
    if (boot == 0) {
      prot.boot = proteins
      t.boot = t.measures
      id.boot = id
      # data for step 2
      surv.boot = survdata
    }
    if (boot > 0) {
      # data for step 1
      keep1 = which(id %in% boot.ids)
      prot.boot = proteins[keep1, ]
      t.boot = t.measures[keep1]
      id.boot = id[keep1]
      # some individuals are sampled more than once...
      # they should be added to the variables and datasets
      # above by changing their identifier
      # get the ids:
      temp1 = as.data.frame(table(boot.ids))
      temp1$boot.ids = as.numeric(as.character(temp1$boot.ids))
      temp1 = temp1[temp1$Freq > 1, ]
      temp1$repetitions = temp1$Freq - 1
      ids.add = rep(temp1$boot.ids, times = temp1$repetitions)
      n.add = length(ids.add)
      u = n + 1
      for (i in ids.add) {
        # NB: extra ids for repeated individuals denoted by u
        # for step 1
        keep2 = which(id == i)
        prot.boot = rbind(prot.boot, proteins[keep2, ])
        t.boot = c(t.boot, t.measures[keep2])
        new.ids = rep(u, length(keep2))
        id.boot = c(id.boot, new.ids)
        u = u+1
      }
      # data for step 2
      temp2 = unique(boot.ids)[order(unique(boot.ids))]
      temp2 = c(temp2, ids.add)
      surv.boot = Surv(time = observed.t[temp2], 
                       event = event[temp2])
    }
    
    # original data for validation:
    # data for step 1
    prot.orig = proteins
    t.orig = t.measures
    id.orig = id
    # data for step 2
    surv.orig = survdata
    
    ################################
    ############ STEP 1 ############
    ################################
    # reconstruct rand int and slopes
    X.boot = matrix(NA, n, 2*n.covs)
    X.orig = matrix(NA, n, 2*n.covs)
    drop.cols = c()
    for (i in 1:n.covs) {
      data.boot = data.frame('y' = prot.boot[ , i],
                             'x' = t.boot,
                             'id' = id.boot)
      data.orig = data.frame('y' = prot.orig[ , i],
                             'x' = t.orig,
                             'id' = id.orig)
      
      # compute predicted random effects
      fixef = as.formula('y ~ x')
      rand.int = 'id'
      rand.slope = 'x'
      temp = get.2ranef(data.boot, data.orig, fixef, rand.int, rand.slope)
      # predicted random effects
      X.boot[ , i] = temp$ranef.train[ , 1] 
      X.orig[ , i] = temp$ranef.test[ , 1]
      if (temp$n.ranef == 2) {
        X.boot[ , n.covs+i] = temp$ranef.train[ , 2]
        X.orig[ , n.covs+i] = temp$ranef.test[ , 2]
      }
      else if (temp$n.ranef == 1) {
        drop.cols = c(drop.cols, n.covs+i)
      }
      rm(temp)
    }
    # drop random slopes when not used
    if (length(drop.cols) > 0) {
      X.boot = X.boot[ , -drop.cols]
      X.orig = X.orig[ , -drop.cols]
    }
      
    ################################
    ############ STEP 2 ############
    ################################
    # fit ridge Cox
    cv.fit = cv.glmnet(x = X.boot, y = surv.boot, family="cox", 
                       standardize = do.stand,
                       alpha = 0, maxit = 1e7)
    
    if (!inherits(cv.fit, 'try-error')) {
      # temporary objects to extract time and status from Surv object:
      temp.boot = as.matrix(surv.boot)
      temp.orig = as.matrix(surv.orig)
      
      ####### C index #######
      # obtain predicted relative risk
      # (for Cox model, type = 'response' gives the predicted relative risk)
      relrisk.boot = predict(cv.fit, newx = X.boot, 
                             s = cv.fit$lambda.min,
                             type = 'response')
      relrisk.orig = predict(cv.fit, newx = X.orig, 
                             s = cv.fit$lambda.min,
                             type = 'response')
      
      # computation of the C index
      # NB: x in concordance.index has to be: a vector of risk predictions
      c.boot = try(concordance.index(x = relrisk.boot, 
                                     surv.time = temp.boot[ , 'time'],
                                     surv.event= temp.boot[ , 'status'], 
                                     method = "noether"))
      c.orig = try(concordance.index(x = relrisk.orig, 
                                     surv.time = temp.orig[ , 'time'],
                                     surv.event= temp.orig[ , 'status'], 
                                     method = "noether"))
      # export C index values
      if (!inherits(c.boot, 'try-error')) out$c.boot[k] = c.boot$c.index
      if (!inherits(c.orig, 'try-error') & boot >0) out$c.orig[k] = c.orig$c.index
      
      ####### time-dependent ROC #######
      # retrieve penalized mle
      pmle.elnet = as.numeric(coef(cv.fit, s = cv.fit$lambda.min))
      linpred.boot = X.boot %*% pmle.elnet
      linpred.orig = X.orig %*% pmle.elnet
      
      for (i in 1:length(t.pred)) {
        ptime = t.pred[i]
        auc.boot = try(survivalROC(Stime = temp.boot[ , 'time'], status = temp.boot[ , 'status'], 
                                   marker = linpred.boot, entry = rep(0, n), 
                                   predict.time = ptime, cut.values = NULL,
                                   method = "NNE", span = 0.25*n^(-0.20)))
        auc.orig = try(survivalROC(Stime = temp.orig[ , 'time'], status = temp.orig[ , 'status'], 
                                   marker = linpred.orig, entry = rep(0, n), 
                                   predict.time = ptime, cut.values = NULL,
                                   method = "NNE", span = 0.25*n^(-0.20)))
        # export time-dependent ROC
        col = 4 + 2*i -1
        if (!inherits(auc.boot, 'try-error') & !is.nan(auc.boot$AUC)) out[k, col] = auc.boot$AUC
        if (!inherits(auc.orig, 'try-error') & !is.nan(auc.orig$AUC) & boot >0) out[k, col+1] = auc.orig$AUC
      }
    }
    # go to next results row
    k = k+1
    
    cat(paste('Boots. sample ', boot, sep = ''))
    cat('\n')
  }
  
  end.time = Sys.time()
    
  resfile = paste(name.part2, a, '_', h, '.RData', sep = '')
  save(out, perc.censoring, 
       start.time, end.time,
       file = resfile)
    
  rm(list = setdiff(ls(), c('a', 'h', 'n.boots', 'n.censoring', 
                            'nsub', 'name.part1', 'name.part2',
                            't.pred', 'get.2ranef', 'do.stand')))
}

