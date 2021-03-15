##################################################################################
## FILE: data_gen_funcs.R
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
library("MASS")

## This data-generating mechanism uses a hierarchical model, where:
## we assume known volumes Vol_j, vol_j for all j (V is total volume of sample, v is volume for qPCR)
## we observe abundances m_j for all j
## and V_{ij} are the qPCR counts for i in 1, 2, ..., q_obs and all j
## and W_{ij} are the br16s relative abundances for all taxa in 1, 2, ..., q and all j
## We assume that:
## (log(p_{1j}/p_{Dj}), ..., log(p_{qj}/p_{Dj})) ~ N(\mu_j, \Sigma_j) for some reference taxa D
## (C_{1j}, ..., C_{qj}) ~ Mult(M_j, p_j), and M_j ~ f(Vol_j)
## (W_{1j}, ..., W_{qj}) ~ Mult(m_j, p_j)
## V_{ij} | p_j ~ Poisson (M_jvol_j/Vol_j p_ij )
## equivalently, we model things on the qPCR scale:
## \beta = mean of the mus (based on data from Hutch collaborators, mean of imputed qpcr)
## \Sigma = coviariance of mus (based on data from Hutch collaborators, covariance of imputed qpcr)
## log \mu_i ~ N_q(\beta, \Sigma) iid for each j, yields \mu_{ij}
## p_{ij} = \mu_{ij}e_i/\sum_i \mu_{ij} e_i
## V_{ij} ~ Poisson(\mu_{ij})
## W_{ij} ~ Mult(M_j, p_{ij})
## ARGS:  beta - hyper parameter, mean for log prob. ratios (q-vector)
##     Sigma - hyper parameter for variance of mus (q-by-q)
##         e - efficiencies (q-vector)
##         q - number of taxa
##         N_subj - number of subjects
##         N_samp - number of samples per subject
##         M - observed br16s reads (N-vector)
##         corr - whether correlation within subjects should be considered
##         epsilon - If corr = T, then the N_subj sized array of sigma_epsilon vals
## RETURNS: X (qPCR), Y (br16s)
data_func <- function(beta, Sigma, e, q, corr = F, sigma_epsilon = rep(1,N_subj),
                      N_subj, N_samp, M) {
  
  ## generate the mus, given the means and covariance;
  ## returns an N x q matrix
  N = N_subj*N_samp
  
  if (corr == F){
    log_mu <- MASS::mvrnorm(N, mu = beta, Sigma = Sigma)
  }
  
  else{
    log_mu <- matrix(NA, nrow = N, ncol = q)
    for (i in 1:N){
      subj <- floor(i/N_samp) + 1
      Sigma_corr <- Sigma + diag(q)*sigma_epsilon[subj]^2
      log_mu[i,] <- MASS::mvrnorm(1, mu = beta, Sigma = Sigma_corr)
    }
  }
  
  mu <- exp(log_mu)
  
  
  ## generate the qPCR
  V <- matrix(NA, nrow = N, ncol = q)
  for (i in 1:q) {
    V[, i] <- rpois(N, mu[, i])
    if (any(is.na(V[, i]))) { # if any are NA, fill with random sample from the others
      non_na <- V[!is.na(V[, i]), i]
      V[is.na(V[, i]), i] <- sample(non_na, sum(is.na(V[, i])), replace = TRUE)
    }
  }
  
  ## generate the br16s reads
  W <- matrix(NA, nrow = N, ncol = q)
  for (j in 1:N) {
    p_j <- e*mu[j, ]/sum(e*mu[j, ])
    W[j, ] <- rmultinom(n = 1, size = M[j], prob = p_j)
  }
  
  ## return
  return(list(V = V, W = W, mu = mu))
}