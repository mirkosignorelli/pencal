###########################
nsub <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(R.version.string)
print(paste('node(s):', Sys.getenv("SLURM_JOB_NODELIST")))
cat('simulation number', nsub, '\n') # simulation number
h = nsub
###########################

on.shark = T

if (on.shark) setwd('/exports/molepi/ms/hdorc3')
if (on.shark == F) {
  rm(list = ls())
  library(msigno)
  clear.screen(); set.cf()
  h = 1
  a = 1
  boot = 1
}

cat('simulation number', nsub, '\n') # simulation number
h = nsub

library(glmnet)
library(survival)
library(survcomp) # C index (concordance.index)
library(survivalROC) # time-dependent AUC
# load function to compute predicted random effects
# on test set using the model estimated on the training set:
source('./functions/mlcmm-BLUP-ranef-v0.3.R')

# standardize predictors in glmnet?
do.stand = T
# step 2: use only u0, u1 or also b0 for predictions?
use.b0s = T

simref = 1
n.cases = 1
name.part1 = paste('data/', simref, '/sim', simref, '_', sep = '')
subf = paste(simref, '.8', sep = '')
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
  markers = markers[-rows.rem, ]
    
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
      prot.boot = markers
      t.boot = t.measures
      id.boot = id
      # data for step 2
      surv.boot = survdata
    }
      
    if (boot > 0) {
      # data for step 1
      keep1 = which(id %in% boot.ids)
      prot.boot = markers[keep1, ]
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
        prot.boot = rbind(prot.boot, markers[keep2, ])
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
    prot.orig = markers
    t.orig = t.measures
    id.orig = id
    # data for step 2
    surv.orig = survdata
      
    ################################
    ############ STEP 1 ############
    ################################
    # reconstruct rand int and slopes
    X.boot = matrix(NA, n, 2*n.prot + n.antib*n.prot)
    X.orig = matrix(NA, n, 2*n.prot + n.antib*n.prot)
    
    col.prot = rep(1:n.prot, each = n.antib)
    for (i in 1:n.prot) {
      cols = which(col.prot == i)
      data.boot = data.frame(prot.boot[ , cols],
                             'x' = t.boot,
                             'id' = id.boot)
      data.orig = data.frame(prot.orig[ , cols],
                             'x' = t.orig,
                             'id' = id.orig)
      y.names = paste('y', 1:n.antib, sep = '')
      names(data.boot)[1:n.antib] = names(data.orig)[1:n.antib] = y.names
      
      # compute predicted random effects
      temp = paste( paste(y.names, collapse = '+'), '~ contrast(x)')
      fixed = as.formula(temp)
      random = ~ 1 + x
      temp = blup.mlcmm(fixed = fixed, random = random,
                 subject = 'id', data = data.boot, newdata = data.orig,
                 maxiter = 1e2)
        
      # predicted random effects
      if (i == 1) {
        X.boot = temp$ranef.orig
        X.orig = temp$ranef.new
      }
      if (i > 1) {
        X.boot = cbind(X.boot, temp$ranef.orig)
        X.orig = cbind(X.orig, temp$ranef.new)
      }
      #print(i)
    }
      
    X.boot = as.matrix(X.boot)
    X.orig = as.matrix(X.orig)
    
    if (!use.b0s) {
      col.keep = which(colnames(X.boot) %in% c('u0', 'u1'))
      X.boot = X.boot[ , col.keep]
      X.orig = X.orig[ , col.keep]
    }
      
    ################################
    ############ STEP 2 ############
    ################################
    # fit elasticnet Cox
    grid = seq(1, 0, by = -0.05)
    n.grid = length(grid)
    fold.ids = rep(1:5, each = n/5)
    models = vector('list', n.grid)
    tuning.matr = matrix(NA, n.grid, 3)
    colnames(tuning.matr) = c('alpha', 'lambda.min', 'deviance')
    
    for (i in 1:n.grid) {
      alpha.el = grid[i]
      temp = try( cv.glmnet(x = X.boot, y = surv.boot, family="cox", 
                            standardize = do.stand,
                            foldid = fold.ids,
                            alpha = alpha.el, maxit = 1e7) )
      tuning.matr[i, ] = c(alpha.el, temp$lambda.min,
                           min(temp$cvm))
      models[[i]] = temp
    }
      
    # optimal combination of (alpha, lambda)
    id.best = which.min(tuning.matr[ ,3])
    cv.fit = models[[id.best]]
    # grid[id.best]
    #plot(coef(cv.fit)[1:200])
    
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
    
    cat(paste('Boots. sample ', boot, ': done', sep = ''))
    cat('\n')
    cat('#####################################')
    cat('\n')
  }
    
  end.time = Sys.time()
  
  resfile = paste(name.part2, a, '_', h, '.RData', sep = '')
  save(out, perc.censoring, 
       start.time, end.time,
       file = resfile)
    
  rm(list = setdiff(ls(), c('a', 'h', 'n.boots', 'n.cases', 'nsub', 
                            'ptime', 'blup.mlcmm', 'compute.sbj.ranef', 
                            'do.stand', 'use.b0s')))
}

