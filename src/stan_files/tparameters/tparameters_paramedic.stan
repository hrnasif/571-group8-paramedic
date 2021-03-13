    //
    vector[q] log_mu;
    // ---- new model -----
    vector[q] p[N_total];
    vector[q] log_mu_t[N_total];
    vector[q_obs] log_mu_tv[N_total];
    
    for (i in 1:N_subj) { 
        if (d > 0) {
            log_mu = beta_0 + (X[i] * beta_1)' + exp(log_Sigma) .* log_mu_tilde[i];
        }
        else {
            log_mu = beta_0 + exp(log_Sigma) .* log_mu_tilde[i];
        }
        // -----new model-----
        for (t in 1:N_samp) {
            log_mu_t[N_samp*(i-1)+t] = log_mu + epsilon_t[N_samp*(i-1)+t];
            if (alpha_sigma > 0 && kappa_sigma > 0) {
                p[N_samp*(i-1)+t] = softmax(log_mu_t[N_samp*(i-1)+t] + log_e);
            }
            else {
                p[N_samp*(i-1)+t] = softmax(log_mu_t[N_samp*(i-1)+t]);
            }
            log_mu_tv[N_samp*(i-1)+t] = head(log_mu_t[N_samp*(i-1)+t], q_obs);
        }
    }