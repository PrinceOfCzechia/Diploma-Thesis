data {
  int<lower=1> N; // number of observations
  int<lower=1> D; // number of ordinal dimensions
  int<lower=2> K; // number of categories per dimension
  int<lower=1> P; // number of predictors
  // Ordinal responses: Y_{i,d} in {1, ..., K}
  // i = 1, ..., N
  // d = 1, ..., D
  array[N, D] int<lower=1, upper=K> Y;
  matrix[N, P] X; // predictors
}

parameters {
  // Regression parameter beta
  // Possibly different for each dimension d=1, ..., D
  matrix[D, P] beta;
  // Thresholds tau_{d,1} < ... < tau_{d,K-1}
  array [D] ordered[K - 1] tau;
  // Cholesky decomposition of D x D latent correlation Sigma
  cholesky_factor_corr[D] L_Sigma;
  // GHK uniform nuisance parameters U_{i,d} ~ U(0,1) iid
  array[N, D] real<lower=0, upper=1> u;
}

model {
  // Priors
  for (d in 1:D) {
    beta[d] ~ normal(0, 5);
    tau[d]  ~ normal(0, 5);
  }
  L_Sigma   ~ lkj_corr_cholesky(4);

  // Likelihood via GHK-style reparameterization
  for (i in 1:N) {
    // Latent standard normal coordinates for observation i
    vector[D] Z_i;

    for (d in 1:D) {
      // Compute conditional mean contribution from previous components Z_i[1, ..., d-1]
      real prev = 0;
      if (d > 1) {
        // GHK formula for contribution from previous components
        // prev = \sum_{j<d} L_Sigma[d,j] * Z_ij
        prev = L_Sigma[d, 1:(d - 1)] * head(Z_i, d - 1);
      }
      // d-th component of the linear predictor for observation i
      real eta_id = row(X, i) * beta[d]';

      // Recursion for dimension d: mu_{i,d} = eta_{i,d} + prev + L_Sigma[d,d] * Z_i[d]
      // We work on the standard-normal scale of Z_i[d]
      real mu_cond  = eta_id + prev;
      real sd_cond = L_Sigma[d, d];

      int k = Y[i, d]; // observed category for (i,d)

      // Standard-normal bounds for z[d] implied by thresholds tau_{k-1}, tau_{k}
      real lower_z_bound;
      real upper_z_bound;

      // tau_0 = -inf, tau_K = +inf
      if (k == 1) {
        lower_z_bound = negative_infinity();
      } else {
        lower_z_bound = (tau[d][k - 1] - mu_cond) / sd_cond;
      }
      if (k == K) {
        upper_z_bound = positive_infinity();
      } else {
        upper_z_bound = (tau[d][k] - mu_cond) / sd_cond;
      }

      // Three cases for truncation of z[d] ~ N(0,1):
      //   1) (-inf, upper_z_bound]     (lowest category)
      //   2) (lower_z_bound, +inf)     (highest category)
      //   3) (lower_z_bound, upper_z_bound]  (interior categories)
      // We transform U_{i,d} ~ U(0,1) to Z_{i,d} with these truncations
      // The Jacobian is log(Phi(upper)-Phi(lower)) for the corresponding interval
      
      real t;
      
      if (k == 1) {
        // Case 1: only upper bound, lower = -inf
        // P(Z_i <= upper_z_bound) = Phi(upper_z_bound) =: b
        // t = b * u  gives t ~ U(0, b)
        // then Z_i = Phi^{-1}(t) is N(0,1) truncated to (-inf, upper_z_bound]
        real b = Phi(upper_z_bound);
        t = b * u[i, d];
        // adding the Jacobian adjustment log(Phi(upper_z_bound))
        target += std_normal_lcdf(upper_z_bound);
      } else if (k == K) {
        // Case 2: only lower bound, upper = +inf
        // P(Z_i <= lower_z_bound) = Phi(lower_z_bound) =: b
        // t = b + (1 - b) * u  gives t ~ U(b, 1)
        // then Z_i = Phi^{-1}(t) is N(0,1) truncated to (lower_z_bound, +inf)
        real b = Phi(lower_z_bound);
        t = b + (1 - b) * u[i, d];
        // adding the Jacobian adjustment log(1 - Phi(lower_z_bound))
        target += std_normal_lccdf(lower_z_bound);
      } else {
        // Case 3: both bounds finite: (lower_z_bound, upper_z_bound]
        // t = Phi(lower_z_bound) + Phi_diff * u  => t ~ U(Phi(lower_z_bound), Phi(upper_z_bound))
        // Z_i = Phi^{-1}(t) is N(0,1) truncated to (lower_z_bound, upper_z_bound]
        real Phi_lower  = Phi(lower_z_bound);
        real Phi_upper  = Phi(upper_z_bound);
        real l_Phi_diff = log_diff_exp(std_normal_lcdf(upper_z_bound), std_normal_lcdf(lower_z_bound));

        t = Phi_lower + exp(l_Phi_diff) * u[i, d];
        // adding the Jacobian adjustment log{Phi(upper_z_bound) - Phi(lower_z_bound)}
        target += l_Phi_diff;
      }
      Z_i[d] = inv_Phi(t);
      // At this point Z_i[d] is a standard normal draw truncated to the bounds
      // implied by Y_{i,d}
    }
  }
}

generated quantities {
  corr_matrix[D] Sigma;
  Sigma = multiply_lower_tri_self_transpose(L_Sigma);
}
