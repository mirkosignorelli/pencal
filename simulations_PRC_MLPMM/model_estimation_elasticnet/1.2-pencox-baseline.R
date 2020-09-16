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
  boot = 1
}

library(glmnet)
library(survival)
library(survcomp) # C index (concordance.index)
library(survivalROC) # time-dependent AUC

# standardize predictors in glmnet?
do.stand = T

simref = 1
n.cases = 1
name.part1 = paste('data/', simref, '/sim', simref, '_', sep = '')
subf = paste(simref, '.2', sep = '')
name.part2 = paste('results/',subf, '/res', subf, '_', sep = '')
a = 1 # no need to repeat a = 1, 2 for this model (same results...)

# number of bootstrap samples
n.boots = 100
# time(s) at which to evaluate tdAUC
t.pred = c(3, 5, 7)

start.time = Sys.time()
  
filename = paste(name.part1, a, '_', h, '.RData', sep = '')
load(filename)
survdata = Surv(time = observed.t, event = event)
# keep only baseline protein values
keep1 = which(t.measures == 0)
id = id[keep1]
markers = markers[keep1, ]
rm(t.measures)
  
out = as.data.frame(matrix(NA, n.boots+1, 2+2+2*length(t.pred)))
names(out) = c('repetition', 'boot', 'c.boot', 'c.orig',
               paste(rep(c('auc.boot.t', 'auc.orig.t'), length(t.pred)), 
                     rep(t.pred, each = 2), sep = ''))

out$repetition = h
k = 1
  
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
    X.boot = markers
    id.boot = id
    surv.boot = survdata
  }
  
  if (boot > 0) {
    X.boot = markers[boot.ids, ]
    id.boot = id[boot.ids]
    surv.boot = survdata[boot.ids]
  }
  
  # original data for validation:
  # data for step 1
  X.orig = markers
  id.orig = id
  # data for step 2
  surv.orig = survdata
  
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
save(out, perc.censoring, start.time, end.time,
       file = resfile)

