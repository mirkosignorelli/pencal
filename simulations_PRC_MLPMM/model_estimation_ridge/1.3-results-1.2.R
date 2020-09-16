rm(list = ls())
library(msigno)
clear.screen(); set.cf()

simref = 1 # 1 to 6
outref = 2 # don't change
subf = paste(simref, '.', outref, sep = '')
name.part2 = paste('results/',subf, '/res', subf, '_', sep = '')
outfile = (paste('results/', simref, '.', outref, '-summaries.RData',
           sep = ''))
n.cases = 1
n.rep = 100

c.list = cens.list = vector('list', 1)
auc3.list = auc5.list = auc7.list = c.list

res = data.frame(matrix(NA, 0, 14))
cens.perc = rep(NA, n.rep)
  
for (h in 1:n.rep) {
  resfile = paste(name.part2, '1_', h, '.RData', sep = '')
  is.ok = try(load(resfile), silent = T)
  if (!inherits(is.ok, 'try-error')) {
    res = rbind(res, out)
    cens.perc[h] = perc.censoring
  }
}

# C index
c.orig = res[res$boot == 0, 1:4]
c.boot = res[res$boot != 0, 1:4]
c.orig$c.naive = c.orig$c.boot
c.boot$optimism = c.boot$c.boot - c.boot$c.orig

# AUC
auc.orig.t3 = res[res$boot == 0, c(1:2, 5:6)]
auc.boot.t3 = res[res$boot != 0, c(1:2, 5:6)]
auc.orig.t5 = res[res$boot == 0, c(1:2, 7:8)]
auc.boot.t5 = res[res$boot != 0, c(1:2, 7:8)]
auc.orig.t7 = res[res$boot == 0, c(1:2, 9:10)]
auc.boot.t7 = res[res$boot != 0, c(1:2, 9:10)]
auc.orig.t3$auc.naive = auc.orig.t3$auc.boot.t3
auc.boot.t3$optimism = auc.boot.t3$auc.boot.t3 - auc.boot.t3$auc.orig.t3
auc.orig.t5$auc.naive = auc.orig.t5$auc.boot.t5
auc.boot.t5$optimism = auc.boot.t5$auc.boot.t5 - auc.boot.t5$auc.orig.t5
auc.orig.t7$auc.naive = auc.orig.t7$auc.boot.t7
auc.boot.t7$optimism = auc.boot.t7$auc.boot.t7 - auc.boot.t7$auc.orig.t7

# compute average optimism correction over bootstrap repetitions
library(dplyr)
mean.opt.c = aggregate(optimism ~ repetition, 
                       data = c.boot,
                       FUN = function(x) mean(x, na.rm = T), drop = F)
mean.opt.auc.t3 = aggregate(optimism ~ repetition, 
                            data = auc.boot.t3,
                            FUN = function(x) mean(x, na.rm = T), drop = F)
mean.opt.auc.t5 = aggregate(optimism ~ repetition, 
                            data = auc.boot.t5,
                            FUN = function(x) mean(x, na.rm = T), drop = F)
mean.opt.auc.t7 = aggregate(optimism ~ repetition, 
                            data = auc.boot.t7,
                            FUN = function(x) mean(x, na.rm = T), drop = F)

# gather naive C, optimism and optimism-corrected C:
c.values = cbind(c.orig[ , c(1, 5)], 'optimism' = mean.opt.c$optimism)
c.values$c.corrected = c.values$c.naive - c.values$optimism
# gather naive AUC, optimism and optimism-corrected AUC:
auc.values.t3 = cbind(auc.orig.t3[ , c(1, 5)], 
                      'optimism' = mean.opt.auc.t3$optimism)
auc.values.t3$auc.corrected = auc.values.t3$auc.naive - auc.values.t3$optimism
auc.values.t5 = cbind(auc.orig.t5[ , c(1, 5)], 
                      'optimism' = mean.opt.auc.t5$optimism)
auc.values.t5$auc.corrected = auc.values.t5$auc.naive - auc.values.t5$optimism
auc.values.t7 = cbind(auc.orig.t7[ , c(1, 5)], 
                      'optimism' = mean.opt.auc.t7$optimism)
auc.values.t7$auc.corrected = auc.values.t7$auc.naive - auc.values.t7$optimism

# save results
c.list[[1]] = as.data.frame(c.values)
auc3.list[[1]] = as.data.frame(auc.values.t3)
auc5.list[[1]] = as.data.frame(auc.values.t5)
auc7.list[[1]] = as.data.frame(auc.values.t7)
cens.list[[1]] = cens.perc

mean.cens = round(sapply(cens.list, function(x) mean(x, na.rm = T), simplify = T), 3)*100
save(mean.cens, c.list,
     auc3.list, auc5.list, auc7.list,
     file = outfile)

