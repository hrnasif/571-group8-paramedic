    // shared prior distributions
    beta_0 ~ normal(0, sigma_beta);
    log_Sigma ~ normal(0, sigma_Sigma + sigma_epsilon);

    if (alpha_sigma > 0 && kappa_sigma > 0) {
        sigma_e ~ inv_gamma(alpha_sigma, kappa_sigma);
        log_e ~ normal(0, sqrt(sigma_e));
    }
    if (d > 0) {
        for (j in 1:q) {
            beta_1[:,j] ~ std_normal();
        }
    }
    if (alpha_phi > 0 && beta_phi > 0) {
        phi ~ gamma(alpha_phi, beta_phi);
    }
