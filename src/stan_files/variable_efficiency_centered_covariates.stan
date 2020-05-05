data{
    int<lower=1> N;
    int<lower=1> q_obs;
    int<lower=1> q;
    int<lower=0> d;
    int<lower=0> V[N,q_obs];
    int<lower=0> W[N,q];
    matrix[N,d] X;
    // hyperparameters
    real sigma_beta;
    real sigma_Sigma;
    real alpha_sigma;
    real kappa_sigma;
}
parameters{
    // first-level parameters
    vector[q] log_mu[N];
    vector[q] log_e;
    // second-level hyperparameters
    vector[q] beta_0;
    matrix[d,q] beta_1;
    vector[q] log_Sigma;
    real<lower=0> sigma_e;
}
transformed parameters{
    vector[q] p[N];
    vector[q_obs] log_mu_v[N];

    for (i in 1:N) {
        p[i] = softmax(log_mu[i] + log_e);
        log_mu_v[i] = head(log_mu[i], q_obs);
    }
}
model {
    // hierarchical model
    beta_0 ~ normal(0, sigma_beta);
    log_Sigma ~ normal(0, sigma_Sigma);

    sigma_e ~ inv_gamma(alpha_sigma, kappa_sigma);
    log_e ~ normal(0, sqrt(sigma_e));

    for (j in 1:q) {
        beta_1[:,j] ~ std_normal();
    }

    for (i in 1:N){
        log_mu[i] ~ normal(beta_0 + (X[i] * beta_1)', exp(log_Sigma));
        V[i] ~ poisson_log(log_mu_v[i]);
        W[i] ~ multinomial(p[i]);
    }
}
generated quantities{
    vector[q] mu[N];
    vector[q] e;
    vector[q] Sigma;
    
    mu = exp(log_mu);
    e = exp(log_e);
    Sigma = exp(log_Sigma);
}
