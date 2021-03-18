setwd("./group8/")
source("./sim-vis/extractinits.R")

# Vary N_samp, q=10, q_obs=7, N_subj=5
q <- 10
qobs <- 7
N_subj <- 5

N_samp <- c(5, 10, 15, 20)

param <- 5
getsim <- function(param){
  folder <- paste("spec-normal-normal-poisson-mult/n_5/n_", 
                  as.character(param),
                  "/q_10/q_obs_7", sep = "")
  
  #truemu <- matrix()
  # N_samp= 5
  data <- readRDS(paste("./sim-data/",folder,"/predict_qpcr_data_jobid_10_ad_0.850000_mt_15_ab_1.rds", sep = ""))
  # Our model
  omod <- readRDS(paste("./sim-data/", folder, "/predict_qpcr_mod_jobid_10_ad_0.850000_mt_15_ab_1.rds", sep = ""))
  osamps <- readRDS(paste("./sim-data/", folder, "/predict_qpcr_samps_jobid_10_ad_0.850000_mt_15_ab_1.rds", sep = ""))
  # Amy's model
  amod <- readRDS(paste("./sim-data/", folder, "/predict_qpcr_mod_jobid_10_ad_0.850000_mt_15_ab_1_Amy.rds", sep = ""))
  asamps <- readRDS(paste("./sim-data/", folder, "/predict_qpcr_samps_jobid_10_ad_0.850000_mt_15_ab_1_Amy.rds", sep = ""))
  
  return(list(data = data, omod = omod, osamps = osamps, amod = amod, asamps = asamps))
}

getrmses <- function(param){
  mods <- getsim(param)
  # True mus
  truemu <- mods$data$mu
  # Model mu estimates
  omu <- mods$omod[1:(N_subj*param*q), 1]
  omu <- matrix(omu, ncol = q, byrow = TRUE)
  amu <- mods$amod[1:(N_subj*param*q), 1]
  amu <- matrix(amu, ncol = q, byrow = TRUE)
  nmu <- inits(data$W, data$V, n_samp = param)
  nmu <- exp(nmu[[1]]$log_mu_tilde)
  
  odiff <- truemu - omu
  ormse <- sqrt(1/(N_samp[1]*N_subj*q)*sum(odiff^2))
  adiff <- truemu - amu
  armse <- sqrt(1/(N_samp[1]*N_subj*q)*sum(adiff^2))
  ndiff <- truemu - nmu
  nrmse <- sqrt(1/(N_samp[1]*N_subj*q)*sum(ndiff^2))
  return(c(ormse, armse, nrmse))
}

rmses <- matrix(nrow = 3, ncol = length(N_samp))
for (i in length(N_samp)){
  param <- N_samp[i]
  rmses[,i] <- getrmses(param)
}

# Credible and prediction intervals
opostsum <- extract_posterior_summaries(omod, osamps, 1:q)
apostsum <- extract_posterior_summaries(amod, asamps, 1:q)


