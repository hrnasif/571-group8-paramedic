data{
    // declares N_subj, N_samp, q_obs, q, d, V, W, X
    // also declares hyperparameters:
    // sigma_beta, sigma_Sigma,
    // alpha_sigma, kappa_sigma,
    // alpha_phi, beta_phi
#include /data/data_paramedic.stan
}
parameters{
    // first-level parameters
    vector[q] log_mu_tilde[N_subj*N_samp];
    // declares shared parameters log_e,
    // beta_0, beta_1, log_Sigma,
    // sigma_e, phi
#include /parameters/parameters_paramedic.stan
}
transformed parameters{
    // declares p, log_mu_tv, log_mu, log_mu_t
#include /tparameters/tparameters_paramedic.stan
}
model {
    // specifies priors on beta_0, beta_1, Sigma,
    // sigma_e, log_e, and phi
#include /model/shared_model_paramedic.stan

    for (i in 1:N_subj) {
        log_mu_tilde[i] ~ std_normal();
        for (t in 1:N_samp) {
            if (alpha_phi > 0 && beta_phi > 0) {
                V[N_samp*(i-1) + t] ~ neg_binomial_2_log(log_mu_tv[N_samp*(i-1) + t], phi[N_samp*(i-1) + t]);
            }
            else {
                V[N_samp*(i-1) + t] ~ poisson_log(log_mu_tv[N_samp*(i-1) + t]);
            }
            W[N_samp*(i-1) + t] ~ multinomial(p[N_samp*(i-1) + t]);
        }
    }
}
generated quantities{
    vector[q] mu_it[N_subj*N_samp];
    vector[q] e;
    vector[q] Sigma;

    mu_it = exp(log_mu_t);
    if (alpha_sigma > 0 && kappa_sigma > 0) {
        e = exp(log_e);
    }
    else {
        e = rep_vector(1, q);
    }
    Sigma = exp(log_Sigma);
}
