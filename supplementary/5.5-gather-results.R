library(foreach)

t.comp = foreach (i = 1:100, .combine = 'rbind') %do% {
  filename = paste('results/5.2/res5.2_1_', i, '.RData', sep = '')
  load(filename)
  time.eval
}

t.comp$time = as.numeric(t.comp$time)
aggregate(time ~ method, FUN = mean, data = t.comp)

results = foreach (i = 1:100, .combine = 'rbind') %do% {
  filename = paste('results/5.4/res5.4_1_', i, '.RData', sep = '')
  load(filename)
  tdauc.vals
}
names(results)[3] = 'PRC LMM'

tdauc2 = subset(results, t.pred == 2)[ , -1]
tdauc3 = subset(results, t.pred == 3)[ , -1]
tdauc4 = subset(results, t.pred == 4)[ , -1]
tdauc5 = subset(results, t.pred == 5)[ , -1]

rounded.mean = function(x) round(mean(x), 3)
aggregate(.~t.pred, data = results, FUN = rounded.mean)

pdf('figs/5.5_tdAUC_validation_set.pdf', width = 11)
par(bty = 'l', mfrow = c(2,2), mar =  c(3,2.5,2,1.5))
boxplot(tdauc2, ylim = c(0.5, 1), main = 'tdAUC, S(2|1)')
boxplot(tdauc3, ylim = c(0.5, 1), main = 'tdAUC, S(3|1)')
boxplot(tdauc4, ylim = c(0.5, 1), main = 'tdAUC, S(4|1)')
boxplot(tdauc5, ylim = c(0.5, 1), main = 'tdAUC, S(5|1)')
dev.off()

pdf('figs/5.5_PRC_vs_joineRML.pdf', width = 10, height = 4)
par(bty = 'l', mfrow = c(1,4), mar =  c(3,2.8,2,1.2), las = 1)
cols = c('darkgoldenrod1', 'cyan2', 'green2')
boxplot(tdauc2, ylim = c(0.5, 1), main = 'tdAUC, t = 2', 
        col = cols, cex.axis = 1.2)
boxplot(tdauc3, ylim = c(0.5, 1), main = 'tdAUC, t = 3', 
        col = cols, cex.axis = 1.2)
boxplot(tdauc4, ylim = c(0.5, 1), main = 'tdAUC, t = 4', 
        col = cols, cex.axis = 1.2)
boxplot(tdauc5, ylim = c(0.5, 1), main = 'tdAUC, t = 5', 
        col = cols, cex.axis = 1.2)
dev.off()

x1 = data.frame(tdauc2$pencal.nopenalty - tdauc2$joineRML,
                tdauc3$pencal.nopenalty - tdauc3$joineRML,
                tdauc4$pencal.nopenalty - tdauc4$joineRML,
                tdauc5$pencal.nopenalty - tdauc5$joineRML)
x2 = data.frame(tdauc2$`PRC LMM` - tdauc2$joineRML,
                tdauc3$`PRC LMM` - tdauc3$joineRML,
                tdauc4$`PRC LMM` - tdauc4$joineRML,
                tdauc5$`PRC LMM` - tdauc5$joineRML)
names(x1) = names(x2) = c('S(2|1)', 'S(3|1)', 'S(4|1)', 'S(5|1)')

pdf('figs/5.5_tdAUC_diff_validation_set.pdf', width = 12)
par(bty = 'l', mfrow = c(1,2), mar =  c(3,4.5,2,1.5))
boxplot(x1, ylab = 'tdAUC: pencal unpenalized - joineRML',
        ylim = c(-0.05, 0.05))
abline(h = 0, col = 'red')
boxplot(x2, ylab = 'tdAUC: pencal ridge - joineRML',
        ylim = c(-0.05, 0.05))
abline(h = 0, col = 'red')
dev.off()
