library(cmdstanr)
library(mvtnorm)

#------------------------------------------------------------------
# Simulate mock data
#------------------------------------------------------------------

set.seed(123)

N = 1000       # observations
D = 3          # dimensions
K = 4          # categories
P = 1          # predictors (single slope)

# Covariates X_{i} (N x P)
X = matrix(rnorm(N * P), nrow = N, ncol = P)

# True regression parameters
beta_true = c(-2.0, 0.0, 1.0)

# True correlation matrix Sigma (D x D)
Sigma_true = matrix(c(
  1.0, 0.7, 0.3,
  0.7, 1.0, 0.5,
  0.3, 0.5, 1.0
), nrow = D, ncol = D, byrow = TRUE)

L_Sigma_true = t(chol(Sigma_true))  # lower Cholesky

# True thresholds tau_1 < tau_2 < tau_3 for 4 categories
tau_true = matrix(
  c(
    -0.5,  0.2,  0.9,   # tau^(1)
    -0.2,  0.4,  1.0,   # tau^(2)
    -0.8,  0.0,  0.8    # tau^(3)
  ),
  nrow = D, ncol = K - 1, byrow = TRUE
)

# Linear predictors
eta = X[, 1, drop = FALSE] %*% t(matrix(beta_true, nrow = D))  # N x D

# Latent Z
Z = mvtnorm::rmvnorm(n = N, mean = rep(0, D), sigma = Sigma_true)
Z = Z + eta

# threshold to ordinal categories 1,...,K
Y = matrix(NA_integer_, nrow = N, ncol = D)
for (d in 1:D) {
  Y[, d] = as.integer(
    cut(
      Z[, d],
      breaks = c(-Inf, tau_true[d, ], Inf),
      labels = FALSE,
      right = TRUE
    )
  )
}

#------------------------------------------------------------------
# Fit the model and inspect posterior
#------------------------------------------------------------------


stan_data = list(
  N = N,
  D = D,
  K = K,
  P = P,
  Y = Y,
  X = X
)

mod = cmdstan_model("ghk_probit_model.stan")

fit = mod$sample(
  data = stan_data,
  seed = 42,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

print(
  fit$summary(
    variables = c(
      "beta",
      "tau",
      "Sigma[1,2]",
      "Sigma[1,3]",
      "Sigma[2,3]"),
    mean, ~quantile(.x, probs = c(0.025, 0.975))
    ),
  n = 15
)

