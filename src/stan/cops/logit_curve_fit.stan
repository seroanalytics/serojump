data {
  int<lower=1> N;          // Number of observations
  vector[N] x;             // Predictor variable (e.g., titre)
  vector[N] y;             // Continuous response variable
}
transformed data {
       
    real x_min = min(x);
    real x_max = max(x);
    real step = (x_max - x_min) / 99;
    real midpoint = (x_min + x_max) / 2;
}

parameters {
  real<lower=0, upper = 1> L;         // Upper asymptote (maximum response)
  real<lower = 0> k;                  // Steepness of the curve
  real x0;                 // Midpoint (inflection point)
  real<lower=0> sigma;     // Standard deviation of errors
}

model {
  vector[N] y_hat;
  
  // Logistic function for curve fitting
  for (n in 1:N) {
    y_hat[n] = L * (1 - 1 / (1 + exp(-k * (x[n] - x0))));
  }

  // Likelihood: Assume normal residuals
  y ~ normal(y_hat, sigma);

  // Priors
  L ~ uniform(0, 1);      // Prior for upper asymptote
  k ~ normal(0, 5);        // Prior for steepness
  x0 ~ normal(midpoint, midpoint / 2);       // Prior for midpoint
  sigma ~ normal(0, 1);    // Prior for standard deviation
}
generated quantities {
    vector[100] y_hat_rel;
    vector[100] y_hat_new;
    vector[100] y_prot;
    vector[100] x_new;
    
    // Posterior prediction
    for (i in 1:100) {
      x_new[i] = x_min + (i - 1) * step;
      y_prot[i] =  1 / (1 + exp(-k * (x_new[i] - x0)));
      y_hat_new[i] = L * (1 - 1 / (1 + exp(-k * (x_new[i] - x0))));
      y_hat_rel[i] = y_hat_new[i] / y_hat_new[1];
    }
}