debug = F
if (debug) {
  library(lme4)
  library(nlme)
  data.train = sleepstudy[1:100, ]
  data.test = sleepstudy[101:180, ]
  fixef = as.formula(Reaction ~ Days)
  rand.int = 'Subject'
  rand.slope = 'Days'
  temp1 = get.1ranef(data.train, data.test, fixef, rand.int)
  temp1
  temp2 = get.2ranef(data.train, data.test, fixef, rand.int, rand.slope)
  temp2
}

get.2ranef = function(data.train, data.test, fixef, rand.int, rand.slope) {
  require(nlme)
  ranef = as.formula(paste('~', rand.slope, '|', rand.int))
  lmm = try( lme(fixed = fixef, random = ranef, data = data.train),
             silent = T)
  if (inherits(lmm, 'try-error')) { # retry with stricter control
    lmm = try( lme(fixed = fixef, random = ranef, data = data.train,
                   control = list(maxIter = 200, msMaxIter = 200,
                                  tolerance = 1e-8, niterEM = 200,
                                  msTol = 1e-8)),
               silent = T)
  }
  
  if (!inherits(lmm, 'try-error')) {
    n.ranef = 2
    Z = model.matrix(as.formula(paste('~', rand.slope)), 
                     data = data.test)
  }
  if (inherits(lmm, 'try-error')) {
    ranef = as.formula(paste('~ 1|', rand.int))
    lmm = try( lme(fixed = fixef, random = ranef, data = data.train) )
    if (inherits(lmm, 'try-error')) { # retry with stricter control
      lmm = try( lme(fixed = fixef, random = ranef, data = data.train,
                     control = list(maxIter = 200, msMaxIter = 200,
                                    tolerance = 1e-8, niterEM = 200,
                                    msTol = 1e-8)) )
    }
    n.ranef = 1
    Z = model.matrix(~ 1, data = data.test)
  } 
  # extract random effects on training set:
  ranef.train = as.matrix(ranef(lmm))
  # compute random effects on test set:
  id.test = data.test[,rand.int]
  n = length(unique(id.test))
  ranef.test = matrix(NA, n, n.ranef)
  D.hat = getVarCov(lmm, type = 'random.effects')
  beta.hat = fixef(lmm)
  sigma2.hat = summary(lmm)$sigma^2
  X = model.matrix(fixef[-2], data = data.test)
  y = data.test[ , as.character(fixef[[2]])]
  for (i in 1:n) {
    id = unique(id.test)[i]
    rows = which(id.test == id)
    Xi = X[rows, , drop = FALSE] #drop prevents conversion to vector when rows has length 1
    yi = y[rows]
    Zi = Z[rows, , drop = FALSE]
    I.matr = diag(1, length(rows), length(rows))
    Vi = Zi %*% D.hat %*% t(Zi) + sigma2.hat * I.matr
    temp = yi - Xi %*% beta.hat
    ranef.test[i, ] = D.hat %*% t(Zi) %*% solve(Vi) %*% temp
  }
  out = list('ranef.train' = ranef.train, 'ranef.test' = ranef.test,
             'lmm.train' = lmm, 'n.ranef' = n.ranef)
  return(out)
}

get.1ranef = function(data.train, data.test, fixef, rand.int) {
  require(nlme)
  ranef = as.formula(paste('~ 1 |', rand.int))
  lmm = try( lme(fixed = fixef, random = ranef, data = data.train),
             silent = T)
  if (inherits(lmm, 'try-error')) { # retry with stricter control
    lmm = try( lme(fixed = fixef, random = ranef, data = data.train,
                   control = list(maxIter = 200, msMaxIter = 200,
                                  tolerance = 1e-8, niterEM = 200,
                                  msTol = 1e-8)),
               silent = T)
  }
  if (!inherits(lmm, 'try-error')) {
    # extract random effects on training set:
    ranef.train = as.matrix(ranef(lmm))
    # compute random effects on test set:
    id.test = data.test[ , rand.int]
    n = length(unique(id.test))
    ranef.test = matrix(NA, n, 1)
    D.hat = getVarCov(lmm, type = 'random.effects')
    beta.hat = fixef(lmm)
    sigma2.hat = summary(lmm)$sigma^2
    X = model.matrix(fixef[-2], data = data.test)
    y = data.test[ , as.character(fixef[[2]])]
    Z = model.matrix(~ 1, data = data.test)
    for (i in 1:n) {
      id = unique(id.test)[i]
      rows = which(id.test == id)
      Xi = X[rows, , drop = FALSE] #drop prevents conversion to vector when rows has length 1
      yi = y[rows]
      Zi = Z[rows, , drop = FALSE]
      I.matr = diag(1, length(rows), length(rows))
      Vi = Zi %*% D.hat %*% t(Zi) + sigma2.hat * I.matr
      temp = yi - Xi %*% beta.hat
      ranef.test[i, ] = D.hat %*% t(Zi) %*% solve(Vi) %*% temp
    }
  }
  if (inherits(lmm, 'try-error')) {
    ranef.train = matrix(0, n, 1)
    ranef.test = matrix(0, n, 1)
  }
  out = list('ranef.train' = ranef.train, 'ranef.test' = ranef.test,
             'lmm.train' = lmm)
  return(out)
}
