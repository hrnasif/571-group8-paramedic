    // ---individual size---
    int<lower=1> N_subj;
    // ---observed times for each individual---
    int<lower=1> N_samp;
    // ---total number of observation---
    int<lower=1> N_total;
    // num taxa with observed absolute abundance
    int<lower=1> q_obs;
    // overall num taxa
    int<lower=1> q;
    // num covariates
    int<lower=0> d;
    // ----observed absolute abundance----
    int<lower=0> V[N_total,q_obs];
    // ----observed rel. abundance-----
    int<lower=0> W[N_total,q];
    // ----feature matrix-----
    matrix[N_total,d] X;
    // hyperparameters that are fixed
    real sigma_beta;
    real sigma_Sigma;
    // 0 for both = efficiency-naive model
    // otherwise, fit varying-efficiency model
    real<lower=0> alpha_sigma;
    real<lower=0> kappa_sigma;
    // 0 for both = Poisson model
    // otherwise, fit negative binomial
    real<lower=0> alpha_phi;
    real<lower=0> beta_phi;
    // hyperparameter for epsilon
    real<lower=0> sigma_epsilon[N_subj];
    