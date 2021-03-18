inits <- function(W, V, X = V[, c(1:2), drop = FALSE], n_samp, k = 0,
                  n_iter = 10, n_burnin = 8000, n_chains = 1, stan_seed = 4747,
                  centered = FALSE, inits_lst = NULL,
                  sigma_beta = sqrt(50), sigma_Sigma = sqrt(50), alpha_sigma = 2, kappa_sigma = 1, 
                  alpha_phi = 0, beta_phi = 0, sigma_xi = 1,
                  ...){
  colnames(W)[1] <- "subject_id"
  colnames(V)[1] <- "subject_id"
  # ---------------------------
  # pre-processing and warnings
  # ---------------------------
  pre_processed_lst <- make_paramedic_tibbles(W, V, X, n_samp, k, inits_lst, sigma_beta, sigma_Sigma, alpha_sigma, kappa_sigma)
  W_mat <- pre_processed_lst$w_mat
  V_mat <- pre_processed_lst$v_mat
  X_mat <- pre_processed_lst$x_mat
  sigma_epsilon <- pre_processed_lst$sigma_epsilon_arry
  # ----------------------------------------
  # set up the data and initial values lists
  # ----------------------------------------
  data_inits_lst <- make_paramedic_stan_data(n_samp, k, W_mat, V_mat, X_mat, inits_lst, 
                                             sigma_beta, sigma_Sigma, 
                                             alpha_sigma, kappa_sigma, 
                                             alpha_phi, beta_phi, sigma_epsilon,
                                             n_chains, centered)
  data_lst <- data_inits_lst$data_lst
  inits_lst <- data_inits_lst$inits_lst
  
  return(inits_lst)
}