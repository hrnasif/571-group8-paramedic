data{
    int<lower=1> N;
    int<lower=1> q_obs;
    int<lower=1> q;
    int<lower=0> V[N,q_obs];
    int<lower=0> W[N,q];
    // hyperparameters
    real sigma_beta;
    real sigma_Sigma;
}
parameters{
    vector[q] mu[N];
    vector[q] beta_0;
    vector[p] beta_1[q];
    vector<lower=0>[q] Sigma;
}
model {
    // hierarchical model
    beta_0 ~ normal(0, sigma_beta);
    Sigma ~ lognormal(0, sigma_Sigma);

    for (j in 1:q) {
        beta_1[j] ~ std_normal();
    }

    for (i in 1:N){
        for (j in 1:q) {
            mu[i,j] ~ lognormal(beta_0[j] + X[i] * beta_1[j], Sigma);
        }
        V[i] ~ poisson(head(mu[i],q_obs));
        W[i] ~ multinomial((mu[i])/sum(mu[i]));
    }
}
