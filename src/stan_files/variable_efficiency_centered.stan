data{
    int<lower=1> N;
    int<lower=1> q_obs;
    int<lower=1> q;
    int<lower=0> V[N,q_obs];
    int<lower=0> W[N,q];
    // hyperparameters
    real sigma_beta;
    real sigma_Sigma;
    real alpha_sigma;
    real kappa_sigma;
}
parameters{
    vector[q] log_mu[N];
    vector[q] beta;
    row_vector<lower=0>[q] Sigma;
    vector<lower=0>[q] e;
    real<lower=0> sigma_e;
}
model {
    // model
    beta ~ normal(0, sigma_beta);
    Sigma ~ lognormal(0, sigma_Sigma);

    sigma_e ~ inv_gamma(alpha_sigma, kappa_sigma);
    e ~ lognormal(0, sqrt(sigma_e));

    for (i in 1:N){
        mu[i] ~ lognormal(beta, Sigma);
        V[i] ~ poisson(head(mu[i],q_obs));
        W[i] ~ multinomial((e .* mu[i])/sum(e .* mu[i]));
    }
}
