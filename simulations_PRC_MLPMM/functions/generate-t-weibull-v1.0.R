# version 1, 24/6/19

generate.t.weibull = function(n, lambda, nu, beta, X, seed = 1) {
  # this formula has been taken from Bender et al. (2005, Statistics 
  # in Medicine), page 1717
  if (lambda <=0 | nu <= 0 | n <=0) stop('lambda, nu and n should be positive!')
  set.seed(seed)
  unif = runif(n, min = 0, max = 1)
  lin.pred = X %*% beta
  if (ncol(lin.pred) > 1) stop('error in computation of the linear predictor')
  temp = - log(unif) / (lambda * exp(lin.pred))
  time = temp^(1/nu)
}

generate.censoring.times = function(n, min, max) {
  out = runif(n, min = min, max = max)
}

moments.weibull = function(lambda, nu) {
  # mean and variance from Bender et al. (2005, Statistics 
  # in Medicine), page 1717
  # median: deduced by setting S(t) = 0.5
  mean = 1/((lambda)^nu)*gamma(1/nu+1)
  median = (-log(0.5)/lambda)^(1/nu)
  var = 1/((lambda^2)^nu) * ( gamma(2/nu+1) - (gamma(1/nu+1))^2 )
  out = list('mean' = mean, 'median' = median, 'std.dev' = sqrt(var))
  return(out)
}