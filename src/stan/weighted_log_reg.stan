data {
    int<lower=0> N;        // Number of unique bins
    vector[N] quantiles;      // Predictor (binned values)
    array[N] int adj_total; // Total count per bin
    array[N] int count_1; // Number of positives per bin
}
parameters {
    real beta_0;   // Intercept
    real beta_1;   // Slope coefficient
 }
transformed parameters {
    vector[N] p = inv_logit(beta_0 + beta_1 * quantiles);
}
model {
  // Priors
  beta_0 ~ normal(0, 5);
  beta_1 ~ normal(0, 5);
  
  // Binomial likelihood (weighted logistic regression)
  count_1 ~ binomial(adj_total, p);
}

generated quantities {
  vector[N] prop_pred;
  
  for (i in 1:N) {
    prop_pred[i] = inv_logit(beta_0 + beta_1 * quantiles[i]); // Predicted proportion
  }

  real k = beta_1;
  real x_0 = -beta_0 / beta_1;
}