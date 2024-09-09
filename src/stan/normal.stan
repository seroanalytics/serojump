functions {
  real custom_p(real esttitreExp, real alpha, real beta, real gamma, real k, real p_base) {
    // Calculate r
    real r = beta * (esttitreExp - alpha);
    
    // Calculate the expression inside the p formula
    real term = (r / pow(1 + abs(r)^k, 1 / k)) * 0.5 + 0.5;
    
    // Calculate p
    real p = gamma * term + p_base;
    
    return p;
  }
}
data {
  int<lower=0> N;          // Number of data points
  vector[N] X;          // Predictor matrix
  vector<lower=0, upper=1>[N] y;  // Outcome variable [0, 1]
  int N_val;
  vector[N_val] X_val;
  real lb_alpha;
}

parameters {
  real<lower = 0> sigma;              // Intercept (for the normal distribution
  
  // 
  real<lower = -10, upper = 10> beta_i;
  real<lower = lb_alpha, upper = 10> alpha;
  real<lower = 0, upper = 1> gamma;
  real<lower = 1, upper = 3> k;
  real<lower = 0, upper = 1> p_base;

}

transformed parameters {
  vector[N] mu;            // Mean of the beta distribution

  for (n in 1:N) {
    mu[n] = custom_p(X[n], alpha, beta_i, gamma, k, p_base);  // Logistic link function
  }
}

model {
  // Priors
  beta_i ~ uniform(-5, 5);           // Prior for coefficients
  alpha ~ uniform(lb_alpha, 10);            // Prior for intercept
  gamma ~ uniform(0, 1);
  k ~ uniform(1, 3);
  p_base ~ uniform(0, 1);
  sigma ~ exponential(1);              // Uniform prior on [0, 1] for precision

  // Likelihood
  y ~ normal(mu, sigma);

}

generated quantities {
  vector[N] y_pred;
  vector[N_val] y_pred_t;
  vector[N_val] y_pred_mean;

  for (n in 1:N) {
    y_pred[n] = normal_rng(mu[n], sigma);
  }


  for (n in 1:N_val) {
    y_pred_mean[n] = custom_p(X_val[n], alpha, beta_i, gamma, k, p_base); 
    y_pred_t[n] = normal_rng(y_pred_mean[n], sigma);
  }


}