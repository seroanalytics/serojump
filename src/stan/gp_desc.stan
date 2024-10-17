functions {
  real monotonic_transform(real x) {
    return -exp(x);  // Enforces the function to be monotonically decreasing
  }
}

data {
  int<lower=0> N;          // Number of data points
  vector[N] x;             // Input variable (predictor)
  vector[N] y;             // Outcome variable
  real<lower=0> length_scale;  // Length scale for the GP
}

parameters {
  vector[N] f_raw;         // Latent function values before applying monotonicity transform
  real<lower=0> sigma;     // Noise standard deviation
}

transformed parameters {
  vector[N] f;             // Latent function values after monotonic transform
  cov_matrix[N] K;         // Covariance matrix

  for (n in 1:N)
    f[n] = monotonic_transform(f_raw[n]);

  // Covariance matrix
  for (i in 1:N) {
    for (j in 1:N) {
      K[i, j] = exp(-0.5 * square((x[i] - x[j]) / length_scale));
    }
  }
  K = K + diag_matrix(rep_vector(square(sigma), N));  // Add noise term to diagonal
}

model {
  // Priors
  f_raw ~ normal(0, 1);   // Latent GP values have a standard normal prior

  // Likelihood
  y ~ multi_normal_cholesky(f, cholesky_decompose(K));  // Multivariate normal distribution
}

generated quantities {
  vector[N] y_pred;          // Posterior predictive distribution
  vector[N] grad_f;          // Gradient of the GP with respect to x

  matrix[N, N] K_inv = inverse(K);  // Inverse of the covariance matrix

  for (n in 1:N) {
    real grad_sum = 0.0;
    for (i in 1:N) {
      grad_sum += K_inv[n, i] * (f_raw[i] / length_scale^2) * (x[n] - x[i]);
    }
    grad_f[n] = grad_sum * monotonic_transform(f_raw[n]);  // Apply monotonic transform to the gradient
  }

  // Generate posterior predictive distribution
  y_pred = multi_normal_rng(f, K);
}