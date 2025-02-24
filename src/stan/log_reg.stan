data {
  int<lower=0> N;          // Number of observations
  int<lower=0> y[N];       // Binary outcome (infection: 1 = Yes, 0 = No)
  vector[N] s;             // Continuous predictor (e.g., titre)
}

parameters {
  real alpha;              // Intercept
  real beta;               // Slope (effect of titre)
}

model {
  // Priors
  alpha ~ normal(0, 5);
  beta ~ normal(0, 2);
  
  // Binomial regression with log link
  y ~ bernoulli_logit(alpha + beta * s);
}

generated quantities {
  vector[N] p_hat;         // Predicted probabilities
  vector[N] rr;            // Relative risk estimates
  
  for (n in 1:N) {
    p_hat[n] = exp(alpha + beta * s[n]);  // Probability of infection at s
  }
  
  // Compute relative risk compared to lowest titre level
  real p_ref = exp(alpha);  // Baseline probability (when s = 0)
  
  for (n in 1:N) {
    rr[n] = p_hat[n] / p_ref;  // Relative risk at s compared to reference
  }
}