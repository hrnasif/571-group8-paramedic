data {
  // ---individual size---
  int<lower=1> N;
  // ---observed times for each individual---
  int<lower=1> T;
  // num taxa with observed absolute abundance
  int<lower=1> q_obs;
  // overall num taxa
  int<lower=1> q;
  // num covariates
  int<lower=0> d;
  // ----observed absolute abundance----
  int<lower=0> V[N*T,q_obs];
  // ----observed rel. abundance-----
  int<lower=0> W[N*T,q];
  // ----feature matrix-----
  matrix[N*T,d] X;
  // hyperparameters that are fixed
  real sigma_beta;
  real sigma_Sigma;
  // 0 for both = efficiency-naive model
  // otherwise, fit varying-efficiency model
  real<lower=0> alpha_sigma;
  real<lower=0> kappa_sigma;
  // ---- longitudinal model-----
  real<lower=0> alpha_epsilon;
  real<lower=0> kappa_epsilon;
  // 0 for both = Poisson model
  // otherwise, fit negative binomial
  real<lower=0> alpha_phi;
  real<lower=0> beta_phi;
}

parameters {
  // shared parameters (no index i, t)
  vector[q] log_e;
  vector[q] epsilon_t;

  // second-level hyperparameters
  // control the mus
    vector[q] beta_0;
    matrix[d,q] beta_1;
    vector[q] log_Sigma;
    // control the es
    real<lower=0> sigma_e;
    // control overdispersion in V
    vector<lower=0>[N] phi;
    //-----new----
    // control the epsilons
    real<lower=0> sigma_epsilon;
}

transformed parameters {
  vector[q] log_mu[N];
  // ---- new model -----
  vector[q] p[N*T];
  vector[q_obs] log_mu_v[N*T];
  vector[q] log_mu_t[N*T];

for (i in 1:N) {
    if (d > 0) {
        log_mu[i] = beta_0 + (X[i] * beta_1)' + exp(log_Sigma) .* log_mu_tilde[i];
    }
    else {
        log_mu[i] = beta_0 + exp(log_Sigma) .* log_mu_tilde[i];
    }
    // -----new model-----
    for (t in 1:T) {
      log_mu_t[T*(i-1)+t] = log_mu[i] + epsilon_t;
      if (alpha_sigma > 0 && kappa_sigma > 0) {
        p[T*(i-1)+t] = softmax(log_mu_t[T*(i-1)+t] + log_e);
        }
      else {
        p[T*(i-1)+t] = softmax(log_mu_t[T*(i-1)+t]);
        }
      log_mu_v[iT*(i-1)+t] = head(log_mu_t[T*(i-1)+t], q_obs);
    }
  }
}

model {
  // shared prior distributions
  beta_0 ~ normal(0, sigma_beta);
  log_Sigma ~ normal(0, sigma_Sigma);

  if (alpha_sigma > 0 && kappa_sigma > 0) {
      sigma_e ~ inv_gamma(alpha_sigma, kappa_sigma);
      log_e ~ normal(0, sqrt(sigma_e));
  }
  // ---- new model-----
  if (alpha_epsilon > 0 && kappa_epsilon > 0) {
    sigma_epsilon ~ inv_gamma(alpha_epsilon, kappa_epsilon);
    epsilon_t ~ normal(0, sqrt(sigma_epsilon))
  }
  //--------------------
  if (d > 0) {
      for (j in 1:q) {
          beta_1[:,j] ~ std_normal();
      }
  }
  if (alpha_phi > 0 && beta_phi > 0) {
      phi ~ gamma(alpha_phi, beta_phi);
  }
}

