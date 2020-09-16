compute.sbj.ranef = function(sbj, data, id.variable,
  name.outcomes, link.mu, link.denom, name.fixefs,
  random, beta, contrs, sigma.u, sigma.b.vec,
  sigma.e.vec) {
  n.outcomes = length(name.outcomes)
  k = which(id.variable == sbj)
  data.sub = data[k, ]
  
  # create Yi
  Yi = c()
  for (i in 1:n.outcomes) {
    temp.Yi = as.numeric(data.sub[name.outcomes[i]][[1]])
    Yi = c(Yi, 
           (temp.Yi - link.mu[i])/link.denom[i] )
  }
  rm(temp.Yi)
  
  # create Xibeta and Zi
  fixef.formula = paste('~', paste(name.fixefs, collapse = '+'))
  Xi<-model.matrix(as.formula(fixef.formula), data.sub)
  Zi.temp<-model.matrix(random, data.sub)
  Zi = Zi.temp
  if (n.outcomes == 1) {
    Xibeta = Xi %*% beta
  }
  if (n.outcomes > 1) {
    beta.temp = beta + c(0, contrs[1])
    Xibeta = Xi %*% beta.temp
    for (i in 2:n.outcomes) {
      beta.temp = beta + c(0, contrs[i])
      Xibeta = rbind(Xibeta, Xi %*% beta.temp)
      # repeat Zi n.outcomes times
      Zi = rbind(Zi, Zi.temp)
    }
  }
  rm(Zi.temp)
  
  # first matrix in the computation
  el1 = t(Zi %*% sigma.u)
  el2 = kronecker(diag(sigma.b.vec), t(rep(1, length(k))))
  mat1 = rbind(el1, el2)
  
  # second matrix
  add1 <-Zi %*% sigma.u %*% t(Zi)
  ntime=dim(Zi)[1]
  rtime=rep(length(k), n.outcomes)
  add2 = diag(rep(sigma.e.vec, rtime))
  library(magic) # for adiag
  add3 = do.call(adiag, mapply(function(a,b)
    matrix(rep(a,b^2),ncol=b),sigma.b.vec,rtime,SIMPLIFY = FALSE))
  sigma.Yi = add1 + add2 + add3
  
  # third matrix
  res = Yi - Xibeta
  
  # computation BLUP
  blup = mat1 %*% solve(sigma.Yi) %*% res
  return(t(blup))
}

blup.mlcmm = function(fixed, random,
                      subject, data, newdata = NULL,
                      maxiter = 1e2) {
  # important: if there is more than 1 covariate in
  # fixed, make sure that fixed = ~ contrast(time) + x2 + x3
  # so that the time variable, with contrasts, is the very first one
  require(lcmm)
  require(foreach)
  n.outcomes = length(all.vars(fixed[[2]]))
  id.variable = data[, subject]
  if (!is.numeric(id.variable)) id.variable = as.numeric(as.factor(id.variable))
  numid = length(unique(id.variable))
  
  # fit model
  mlcmm.fit = multlcmm(fixed = fixed, 
                      random = random, 
                      randomY = T,
                      subject = subject, 
                      data = data, 
                      maxiter = maxiter)
  if (mlcmm.fit$conv == 2) warning('maximum # of iterations reached')
  else if (mlcmm.fit$conv > 2) warning('probable convergence issue')
  
  fit.names = names(mlcmm.fit$best)
  name.outcomes = all.vars(fixed[[2]])
  name.fixefs = all.vars(fixed[-2])
  
  # get positions
  pos.fixefs = which(fit.names %in% name.fixefs)
  pos.contrasts = which(substr(fit.names, 1, 8) == 'contrast')
  pos.sigma.u = which(substr(fit.names, 1, 6) == 'varcov')
  pos.sigma.e = which(substr(fit.names, 1, 7) == 'std.err')
  pos.sigma.b = which(substr(fit.names, 1, 11) == 'std.randomY')
  pos.link.mu = which(substr(fit.names, 1, 8) == 'Linear 1')
  pos.link.denom = which(substr(fit.names, 1, 8) == 'Linear 2')
  
  # retrieve variance matrices
  if (length(pos.sigma.u) != 2) {
    stop('code currently supports only 1 shared random intercept and slope')
  }
  if (length(pos.sigma.u) == 2) {
    sigma.u = matrix(c(1, mlcmm.fit$best['varcov 1'],
                       mlcmm.fit$best['varcov 1'], 
                       mlcmm.fit$best['varcov 2']), 2, 2)
  }
  sigma.e.vec = mlcmm.fit$best[pos.sigma.e]^2
  sigma.b.vec = mlcmm.fit$best[pos.sigma.b]^2
  
  # get betas and links
  beta = c(0, mlcmm.fit$best[pos.fixefs])
  contrs = mlcmm.fit$best[pos.contrasts]
  last.contr = -sum(contrs)
  contrs = c(contrs, last.contr) # add the last element
  link.mu = mlcmm.fit$best[pos.link.mu]
  link.denom = mlcmm.fit$best[pos.link.denom]
  
  # compute estimates of random effects on original data
  ranef.orig = matrix(NA, nrow = numid, ncol = 2+n.outcomes)
  unique.sbjid = unique(id.variable)
  
  ranef.orig = foreach(sbj = unique.sbjid, .combine = rbind) %do%
    compute.sbj.ranef(sbj = sbj, data = data, id.variable = id.variable,
                      name.outcomes = name.outcomes, link.mu = link.mu, 
                      link.denom = link.denom, name.fixefs = name.fixefs,
                      random = random, beta = beta, contrs = contrs, 
                      sigma.u = sigma.u, sigma.b.vec = sigma.b.vec, 
                      sigma.e.vec = sigma.e.vec)
  
  colnames(ranef.orig) = c('u0', 'u1',
       paste('b', 1:n.outcomes, sep = ''))
  rownames(ranef.orig) = unique.sbjid
  ranef.orig = as.data.frame(ranef.orig)
  
  # if newdata provided, compute them also on the new dataset
  if (!is.null(newdata)) {
    if (!identical(names(newdata), names(data))) {
      stop('variable names in newdata are different from
           variable names in data')
    }
    id.variable = newdata[, subject]
    if (!is.numeric(id.variable)) id.variable = as.numeric(as.factor(id.variable))
    numid = length(unique(id.variable))
    ranef.new = matrix(NA, nrow = numid, ncol = 2+n.outcomes)
    unique.sbjid = unique(id.variable)
    
    ranef.new = foreach(sbj = unique.sbjid, .combine = rbind) %do%
      compute.sbj.ranef(sbj = sbj, data = newdata, id.variable = id.variable,
                        name.outcomes = name.outcomes, link.mu = link.mu, 
                        link.denom = link.denom, name.fixefs = name.fixefs,
                        random = random, beta = beta, contrs = contrs, 
                        sigma.u = sigma.u, sigma.b.vec = sigma.b.vec, 
                        sigma.e.vec = sigma.e.vec)
    colnames(ranef.new) = c('u0', 'u1',
                            paste('b', 1:n.outcomes, sep = ''))
    rownames(ranef.new) = unique.sbjid
    ranef.new = as.data.frame(ranef.new)
  }
  # define exports
  out = list('ranef.orig' = ranef.orig)
  if (!is.null(newdata)) out[['ranef.new']] = ranef.new
  out[['mle.mlcmm']] = mlcmm.fit$best
  return(out)
}
