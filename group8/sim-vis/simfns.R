getsim <- function(N_samp, q, qobs){
  folder <- paste("spec-normal-normal-poisson-mult/n_5/n_", 
                  as.character(N_samp),
                  "/q_", as.character(q), "/q_obs_", as.character(qobs), sep = "")
  
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

getsim_misspec <- function(simname){
  folder <- paste(simname, "/n_5/n_10/q_10/q_obs_7", sep = "")
  
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

getsim_misspec_ours <- function(simname){
  folder <- paste(simname, "/n_5/n_10/q_10/q_obs_7", sep = "")
  
  #truemu <- matrix()
  # N_samp= 5
  data <- readRDS(paste("./sim-data/",folder,"/predict_qpcr_data_jobid_10_ad_0.850000_mt_15_ab_1.rds", sep = ""))
  # Our model
  omod <- readRDS(paste("./sim-data/", folder, "/predict_qpcr_mod_jobid_10_ad_0.850000_mt_15_ab_1.rds", sep = ""))
  osamps <- readRDS(paste("./sim-data/", folder, "/predict_qpcr_samps_jobid_10_ad_0.850000_mt_15_ab_1.rds", sep = ""))
  
  return(list(data = data, omod = omod, osamps = osamps))
}

mus <- function(mods, N_samp, q){
  omu <- mods$omod[1:(N_subj*N_samp*q), 1]
  omu <- matrix(omu, ncol = q, byrow = TRUE)
  amu <- mods$amod[1:(N_subj*N_samp*q), 1]
  amu <- matrix(amu, ncol = q, byrow = TRUE)
  nmu <- inits(mods$data$W, mods$data$V, n_samp = N_samp)
  nmu <- exp(nmu[[1]]$log_mu_tilde)
  return(list(omu = omu, amu = amu, nmu = nmu))
}

getrmses <- function(N_samp, q, qobs, mods){
  # True mus
  truemu <- mods$data$mu
  # Model mu estimates
  mus <- mus(mods, N_samp, q)
  
  odiff <- truemu - mus$omu
  ormse <- sqrt(1/(N_subj*N_samp*q)*sum(odiff^2))
  adiff <- truemu - mus$amu
  armse <- sqrt(1/(N_subj*N_samp*q)*sum(adiff^2))
  ndiff <- truemu - mus$nmu
  nrmse <- sqrt(1/(N_subj*N_samp*q)*sum(ndiff^2))
  return(c(ormse, armse, nrmse))
}

getrmspes <- function(N_samp, q, qobs, mods){
  data <- mods$data
  vstar <- data$Vstar[,(qobs + 1):q]
  
  mus <- mus(mods, N_samp, q)
  ov <- mus$omu[,(qobs + 1):q]
  av <- mus$amu[,(qobs + 1):q]
  nv <- mus$nmu[,(qobs + 1):q]
  
  odiff <- vstar - ov
  ormspe <- sqrt(1/(N_subj*N_samp*q)*sum(odiff^2))
  adiff <- vstar - av
  armspe <- sqrt(1/(N_subj*N_samp*q)*sum(adiff^2))
  ndiff <- vstar - nv
  nrmspe <- sqrt(1/(N_subj*N_samp*q)*sum(ndiff^2))
  return(c(ormspe, armspe, nrmspe))
}

getmucoverage <- function(N_samp, q, qobs, mods){
  # True mus
  truemu <- mods$data$mu
  
  # Credible and prediction intervals
  opostsum <- extract_posterior_summaries(mods$omod, mods$osamps, 1:q)
  apostsum <- extract_posterior_summaries(mods$amod, mods$asamps, 1:q)
  
  ocred <- opostsum$cred_intervals
  acred <- apostsum$cred_intervals
  
  sums <- c()
  for(i in 1:q){
    inti <- ocred[,,i]
    truei <- truemu[,i]
    indicatori <- (inti[,1] <= truei) & (inti[,2] >= truei)
    sums <- append(sums, sum(indicatori))
  }
  ocovmu <- 1/(N_subj*N_subj*q)*sum(sums)
  
  sums <- c()
  for(i in 1:q){
    inti <- acred[,,i]
    truei <- truemu[,i]
    indicatori <- (inti[,1] <= truei) & (inti[,2] >= truei)
    sums <- append(sums, sum(indicatori))
  }
  acovmu <- 1/(N_subj*N_subj*q)*sum(sums)
  return(c(ocovmu, acovmu))
}

getvcoverage <- function(N_samp, q, qobs, mods){
  data <- mods$data
  vstar <- data$Vstar[,(qobs + 1):q]
  
  opostsum <- extract_posterior_summaries(mods$omod, mods$osamps, 1:q)
  apostsum <- extract_posterior_summaries(mods$amod, mods$asamps, 1:q)
  
  opred <- opostsum$pred_intervals
  apred <- apostsum$pred_intervals
  
  sums <- c()
  for(i in (qobs+1):q){
    inti <- opred[,,i]
    truei <- vstar[,(i-qobs)]
    indicatori <- (inti[,1] <= truei) & (inti[,2] >= truei)
    sums <- append(sums, sum(indicatori))
  }
  ocovv <- 1/(N_subj*N_samp*(q-qobs))*sum(sums)
  
  sums <- c()
  for(i in (qobs+1):q){
    inti <- apred[,,i]
    truei <- vstar[,(i-qobs)]
    indicatori <- (inti[,1] <= truei) & (inti[,2] >= truei)
    sums <- append(sums, sum(indicatori))
  }
  acovv <- 1/(N_subj*N_samp*(q-qobs))*sum(sums)
  return(c(ocovv, acovv))
}
