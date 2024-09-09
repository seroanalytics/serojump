data {
  int<lower=0> N;          // Number of data points
  vector[N] X;          // Predictor matrix
  vector<lower=0, upper=1>[N] y;  // Outcome variable [0, 1]
}

parameters {
  real beta;          // Coefficients for predictors (for beta part)
  real alpha;              // Intercept (for beta part)
  real<lower=0, upper=1> phi;  // Precision parameter for the beta distribution
  real<lower=0, upper=1> zeta_0; // Probability of y = 0
  real<lower=0, upper=1> zeta_1; // Probability of y = 1
}

transformed parameters {
  vector[N] mu;            // Mean of the beta distribution

  for (n in 1:N) {
    mu[n] = inv_logit(alpha + X[n] * beta);  // Logistic link function
  }
}

model {
  // Priors
  beta ~ normal(0, 5);           // Prior for coefficients
  alpha ~ normal(0, 5);          // Prior for intercept
  phi ~ beta(1, 1);              // Uniform prior on [0, 1] for precision
  zeta_0 ~ beta(2, 2);           // Prior for probability of zero inflation
  zeta_1 ~ beta(2, 2);           // Prior for probability of one inflation

  // Likelihood
  for (n in 1:N) {
    if (y[n] == 0) {
      target += log(zeta_0);  // Contribution for y = 0
    } else if (y[n] == 1) {
      target += log(zeta_1);  // Contribution for y = 1
    } else {
      target += log1m(zeta_0 + zeta_1) +
                beta_lpdf(y[n] | mu[n] * phi, (1 - mu[n]) * phi);  // Beta part
    }
  }
}

generated quantities {
  vector[N] y_pred;

  for (n in 1:N) {
    real p_0 = zeta_0;
    real p_1 = zeta_1;
    real p_beta = 1 - zeta_0 - zeta_1;

    if (bernoulli_rng(p_0)) {
      y_pred[n] = 0;
    } else if (bernoulli_rng(p_1 / (p_1 + p_beta))) {
      y_pred[n] = 1;
    } else {
      y_pred[n] = beta_rng(mu[n] * phi, (1 - mu[n]) * phi);
    }
  }

    for (n in 1:10) {
    real p_0 = zeta_0;
    real p_1 = zeta_1;
    real p_beta = 1 - zeta_0 - zeta_1;

    if (bernoulli_rng(p_0)) {
      y_pred[n] = 0;
    } else if (bernoulli_rng(p_1 / (p_1 + p_beta))) {
      y_pred[n] = 1;
    } else {
      y_pred[n] = beta_rng(mu[n] * phi, (1 - mu[n]) * phi);
    }
  }




}