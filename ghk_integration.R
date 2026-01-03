library(mvtnorm)

ghk_pmvnorm = function(mu, Sigma, a, b, K = 1e4) {
  # mu: D-vector
  # Sigma: D x D covariance, PD matrix
  # a, b: D-vectors of bounds, possibly  +-Inf
  # K: number of GHK draws
  
  D  = length(mu)
  L = t(chol(Sigma)) # Cholesky decomposition lower triangle
  
  h_draws = numeric(K)
  
  for (k in seq_len(K)) {
    x_tilde = numeric(D)
    h = 1.0
    
    for (j in seq_len(D)) {
      # conditional mean mu_{j|<j}(x_tilde_{<j})
      if (j == 1) {
        mu_cond = mu[1]
        sd_cond = L[1, 1]
      } else {
        # contribution from previous components via Cholesky
        mu_cond = mu[j] + sum(L[j, 1:(j - 1)] * x_tilde[1:(j - 1)])
        sd_cond = L[j, j]
      }
      
      lo_std = (a[j] - mu_cond) / sd_cond
      hi_std = (b[j] - mu_cond) / sd_cond
      
      Phi_lo   = pnorm(lo_std)
      Phi_hi   = pnorm(hi_std)
      Phi_diff = Phi_hi - Phi_lo
      
      u = runif(1)
      t = u * Phi_hi + (1 - u) * Phi_lo
      x_tilde[j] = qnorm(t)
      
      h = h * Phi_diff
    }
    h_draws[k] = h
  }
  mean(h_draws)
}


rho = -0.6
Sigma = matrix(c(1, 0.8, 0.6, 0.4,
                 0.8, 1, 0.8, 0.6,
                 0.6, 0.8, 1, 0.8,
                 0.4, 0.6, 0.8, 1), 4, 4)
mu = c(0, 0, 0, 0)
a =  c(-5, -2, 0, -20)
b = c(0, 2, 3, 20)

# GHK estimate
set.seed(42)
ghk_est = ghk_pmvnorm(mu=mu, Sigma=Sigma, a=a, b=b, K=1e4)

# mvtnorm integral
mvtnorm_est = pmvnorm(lower=a, upper=b, mean=mu, sigma=Sigma)

c(ghk = ghk_est, mvtnorm = mvtnorm_est)
