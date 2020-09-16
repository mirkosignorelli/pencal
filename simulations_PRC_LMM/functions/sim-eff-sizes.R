sim.eff.sizes = function(n, abs.range) {
  # simulates from a U(a,b) and adds a sign (p = 0.5)
  abs.val = runif(n, min = abs.range[1], max = abs.range[2])
  signs = 2*rbinom(n, 1, prob = 0.5) - 1
  out = signs*abs.val
  return(out)
}