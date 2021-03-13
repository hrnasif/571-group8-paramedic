  // sample size
  int <lower = 1> N;
  // number of samples
  int <lower=1> N_samples;
  // num taxa with observed absolute abundance
  int<lower=1> q_obs;
  // overall num taxa
  int<lower=1> q;
  // read counts
  int<lower=1> M[N];
  // num covariates
  int<lower=0> d;
  // feature matrix (test data)
  matrix[N,d] X;
  // posterior distributions on sigma_e
  vector[N_samples] sigma_e;
  // posterior distributions on beta_0, beta_1, Sigma, phi
  vector[q] beta_0[N_samples];
  matrix[d, q] beta_1[N_samples];
  vector[q] Sigma[N_samples];
  vector[N] phi[N_samples];
  // 0 for both = efficiency-naive model
  // otherwise, fit varying-efficiency model
  real<lower=0> alpha_sigma;
  real<lower=0> kappa_sigma;
  // 0 for both = Poisson model
  // otherwise, fit negative binomial
  real<lower=0> alpha_phi;
  real<lower=0> beta_phi;

  /* Sylvia's edit
  // number of subjects
  int<lower=1> N_subj;
  // number of samples for each subject
  int<lower=1> N_samp;
  //--- for new model, we need Sigma_epsilon---
  real<lower=0> sigma_epsilon[N_subj];
  */