functions {
  // Custom log likelihood function for logit-normal distribution
  real logit_normal_lpdf(vector y, vector mu, real sigma) {
    return normal_lpdf(logit(y) | mu, sigma) - sum(log(y) + log1m(y));
  }


}

data {
  int<lower=0> N;          // Number of data points
  vector[N] X;          // Predictor matrix
  vector<lower=0, upper=1>[N] y;  // Outcome variable (0 < y < 1)
}

parameters {
  real beta_i;          // Coefficients for predictors
  real alpha;              // Intercept
  real<lower=0> sigma;     // Standard deviation of the error term
}

transformed parameters {
  vector[N] mu;            // Linear predictor for the mean of the normal distribution

  for (n in 1:N) {
    mu[n] = alpha + X[n] * beta_i;  // Linear predictor
  }
}

model {
  // Priors
  beta_i ~ normal(0, 5);       // Prior for coefficients
  alpha ~ normal(0, 5);      // Prior for intercept
  sigma ~ cauchy(0, 2.5);    // Prior for standard deviation

  // Likelihood using custom logit-normal log likelihood
  y ~ logit_normal(mu, sigma);
}

generated quantities {
  vector[N] y_pred;          // Posterior predictive distribution

  for (n in 1:N) {
    y_pred[n] = inv_logit(normal_rng(mu[n], sigma));  // Generate from posterior
  }
}