data {
  int<lower=1> N;               // Number of observations (10 points)
  vector[N] x;                  // Predictor variable (e.g., titre)
  vector[N] y;                  // Continuous response variable
  vector<lower=0>[N] weights;   // Weights for each observation
}

parameters {
  real<lower=0, upper=1> L;      // Upper asymptote (max response)
  real k;                        // Steepness of the curve
  real x0;                       // Midpoint (inflection point)
  real<lower=0> sigma;           // Standard deviation of errors
}

transformed parameters {
  vector[N] y_hat;
  
  // Compute the logistic function for each observation
  for (n in 1:N) {
    y_hat[n] = L * (1 - 1 / (1 + exp(-k * (x[n] - x0))));
  }
}
model {
  // Priors
  L ~ uniform(0, 1);            // Upper asymptote between 0 and 1
  k ~ normal(0, 5);             // Steepness prior
  x0 ~ normal(0, 5);            // Midpoint prior
  sigma ~ normal(0, 1);         // Standard deviation prior
  
  // Weighted likelihood
  for (n in 1:N) {
    target += weights[n] * student_t_lpdf(y[n] | 3, y_hat[n], sigma);
  }
}
generated quantities {
  vector[N] y_hat_new;
  
  for (n in 1:N) {
    y_hat_new[n] = 1 / (1 + exp(-k * (n - x0)));
  }
}
