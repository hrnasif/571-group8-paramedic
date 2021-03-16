##################################################################################
## FILE: data_generator.R
##
## CREATED: 11 Jan 2018 by Brian Williamson
## EDITED: 15 MAR 2021 by Hassan Nasif
##
## PURPOSE: Generate data based on simulation setting
##
## INPUTS: 
##         
## OUTPUTS: a stan model file saved as a .rds
##################################################################################

# create data
data_generator <- function(N_subj, N_samp, num_taxa, num_qpcr, seed, 
                           hyper_mean_mu, hyper_cov_mu,
                           hyper_sigma, corr_within = F, hyper_sigma_epsilon,
                           hyper_m_min, hyper_m_max, use_most_abundant) {
  
  sample_size = N_subj * N_samp
  
  # get the efficiencies
  e <- exp(rnorm(num_taxa, mean = 0, sd = hyper_sigma))
  
  # get the copy number
  # cpy <- rep(1, sample_size)
  # get the m's (observed number of reads for each subject)
  m <- sample(hyper_m_min:hyper_m_max, sample_size, replace = TRUE)
  
  # set up a return list
  ret = list(Vstar = matrix(NA, sample_size, num_taxa), W = matrix(NA, sample_size, num_taxa),
             V = matrix(NA, sample_size, num_qpcr), N = as.integer(sample_size), q = as.integer(num_taxa),
             q_obs = as.integer(num_qpcr), 
             mu = matrix(NA, sample_size, num_taxa), e = e, sigma = hyper_sigma, beta = hyper_mean_mu, 
             Sigma = hyper_cov_mu, sigma_epsilon = hyper_sigma_epsilon,
             a = hyper_m_min, b = hyper_m_max, m = m)
  
  # set the seed
  set.seed(seed)
  
  # generate the data
  tmp <- data_func(hyper_mean_mu, hyper_cov_mu, e, num_taxa, corr_within, hyper_sigma_epsilon,
                   N_subj, N_samp, m)
  ret$mu <- tmp$mu
  # knock out some qPCR information; always make it the final column(s) of Xstar, mu, Y
  ret$Vstar <- tmp$V
  ret$W <- tmp$W
  if (use_most_abundant) {
    ## use only the q^obs most abundant taxa in algorithm (scramble the columns to match)
    taxa_avg_abundance <- colMeans(tmp$W/m)
    ordered_avg_abundance <- order(taxa_avg_abundance, decreasing = TRUE)
    ret$W <- tmp$W[, ordered_avg_abundance]
    ret$Vstar <- tmp$V[, ordered_avg_abundance]
    ret$beta <- beta[ordered_avg_abundance]
    ret$e <- e[ordered_avg_abundance]
    ret$mu <- tmp$mu[, ordered_avg_abundance]
  }
  if (num_qpcr < num_taxa) {
    for (i in 1:num_taxa) {
      if (i <= num_qpcr) {
        ret$V[, i] <- ret$Vstar[, i]
      }
    }
  } else {
    ret$V <- ret$Vstar
  }
  subject_ids <- rep(1:N_subj, each=N_samp)
  collection_times <- rep(1:N_samp, times=N_subj)
  ret$W <- cbind(subject_ids, collection_times, ret$W)
  ret$V <- cbind(subject_ids, collection_times, ret$V)
  
  return(ret)
}