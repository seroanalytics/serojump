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
  real ub_alpha;
  // heir
  int N_g1;
  array[N] int X_g1;
  array[N_g1] real prop_g1;
}
transformed data{
  real k = 1;
  real p_base = 0;
}
parameters {
  real<lower = 0> sigma;              // Intercept (for the normal distribution
  real<lower = -10, upper = 10> beta_i;
  real<lower = -10, upper = 10> alpha;
  real<lower = -10, upper = 10> gamma;

  // heirarchicial
  array[N_g1] real z_alpha;
  array[N_g1] real z_beta;
  array[N_g1] real z_gamma;

  real<lower = 0> alpha_sigma;
  real<lower = 0> beta_sigma;
  real<lower = 0> gamma_sigma;

}

transformed parameters {
  vector[N] mu;            // Mean of the beta distribution

  for (n in 1:N) {
    real alpha_reg = inv_logit(alpha + z_alpha[X_g1[n]] * alpha_sigma) * (ub_alpha - lb_alpha) + (lb_alpha);
    real beta_i_ref = inv_logit(beta_i + z_beta[X_g1[n]] * beta_sigma) * -2;
    real gamma_ref = inv_logit(gamma + z_gamma[X_g1[n]] * gamma_sigma);

    mu[n] = custom_p(X[n], alpha_reg, beta_i_ref, gamma_ref, k, p_base);  // Logistic link function
  }
}

model {
  // Priors
  beta_i ~ normal(0, 1.68);           // Prior for coefficients
  alpha ~ normal(0, 1.68);          // Prior for intercept
  gamma ~ normal(0, 1.68);
  //p_base ~ uniform(0, 1);
  sigma ~ exponential(10);              // Uniform prior on [0, 1] for precision

  // Likelihood
  y ~ normal(mu, sigma);

  z_alpha ~ std_normal();
  z_beta ~ std_normal();
  z_gamma ~ std_normal();

  alpha_sigma ~ exponential(3);
  beta_sigma ~ exponential(3);
  gamma_sigma ~ exponential(3);

}

generated quantities {
  vector[N] y_pred;
  matrix[N_val, N_g1] y_pred_t;
  matrix[N_val, N_g1] y_pred_mean_cov;
  vector[N_val] y_pred_mean;

  for (n in 1:N) {
    y_pred[n] = normal_rng(mu[n], sigma);
  }

  real alpha_reg_mean_out = 0;
  real beta_i_ref_mean_out = 0;
  real gamma_ref_mean_out = 0;

  for (j in 1:N_g1) {
    real alpha_reg_out = inv_logit(alpha + z_alpha[j] * alpha_sigma) * (ub_alpha - lb_alpha) + (lb_alpha);
    real beta_i_ref_out = inv_logit(beta_i + z_beta[j] * beta_sigma) * -2;
    real gamma_ref_out = inv_logit(gamma + z_gamma[j] * gamma_sigma);

    alpha_reg_mean_out += alpha_reg_out * prop_g1[j];
    beta_i_ref_mean_out += beta_i_ref_out * prop_g1[j];
    gamma_ref_mean_out += gamma_ref_out * prop_g1[j];
}

  for (n in 1:N_val) {
    real alpha_reg_mean = 0;
    real beta_i_ref_mean = 0;
    real gamma_ref_mean = 0;
    
    for (j in 1:N_g1) {
      real alpha_reg = inv_logit(alpha + z_alpha[j] * alpha_sigma) * (ub_alpha - lb_alpha) + (lb_alpha);
      real beta_i_ref = inv_logit(beta_i + z_beta[j] * beta_sigma) * -2;
      real gamma_ref = inv_logit(gamma + z_gamma[j] * gamma_sigma);

      alpha_reg_mean += alpha_reg * prop_g1[j];
      beta_i_ref_mean += beta_i_ref * prop_g1[j];
      gamma_ref_mean += gamma_ref * prop_g1[j];

      y_pred_mean_cov[n, j] = custom_p(X_val[n], alpha_reg, beta_i_ref, gamma_ref, k, p_base); 
      y_pred_t[n, j] = normal_rng(y_pred_mean_cov[n, j], sigma);
    }
    y_pred_mean[n]  = custom_p(X_val[n], alpha_reg_mean, beta_i_ref_mean, gamma_ref_mean, k, p_base); 
  }
}