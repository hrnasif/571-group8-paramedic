#source("./sim-vis/simfns.R")
source("./sim-vis/extractinits.R")

simnames <- c("misspec-gamma-gamma-poisson-mult",
"misspec-halft-gamma-poisson-mult",
"misspec-gamma-halft-poisson-mult",
"misspec-halft-halft-poisson-mult",
"misspec-normal-gamma-poisson-mult",
"misspec-normal-halft-poisson-mult",
"misspec-gamma-normal-poisson-mult",
"misspec-halft-normal-poisson-mult",
"misspec-normal-normal-poisson-dirimult",
"misspec-normal-normal-negbin-mult",
"misspec-normal-normal-negbin-dirimult",
"spec-normal-normal-poisson-mult"
)
simnames_short <- c("gamma, gamma, poisson, mult",
                    "halft, gamma, poisson, mult",
                    "gamma, halft, poisson, mult",
                    "halft, halft, poisson, mult",
                    "normal, gamma, poisson, mult",
                    "normal, halft, poisson, mult",
                    "gamma, normal, poisson, mult",
                    "halft, normal, poisson, mult",
                    "normal, normal, poisson, dirimult",
                    "normal, normal, negbin, mult",
                    "normal, normal, negbin, dirimult",
                    "normal, normal, poisson, mult"
)
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
  return(list(omu=omu))
}

getrmses <- function(N_samp, q, qobs, mods){
  # True mus
  truemu <- mods$data$mu
  # Model mu estimates
  mus <- mus(mods, N_samp, q)
  
  odiff <- truemu - mus$omu
  ormse <- sqrt(1/(N_subj*N_samp*q)*sum(odiff^2))
  return(c(ormse))
}

getrmspes <- function(N_samp, q, qobs, mods){
  data <- mods$data
  vstar <- data$Vstar[,(qobs + 1):q]
  
  mus <- mus(mods, N_samp, q)
  ov <- mus$omu[,(qobs + 1):q]
  
  odiff <- vstar - ov
  ormspe <- sqrt(1/(N_subj*N_samp*q)*sum(odiff^2))
  
  return(c(ormspe))
}

getmucoverage <- function(N_samp, q, qobs, mods){
  # True mus
  truemu <- mods$data$mu
  
  # Credible and prediction intervals
  opostsum <- extract_posterior_summaries(mods$omod, mods$osamps, 1:q)
  
  ocred <- opostsum$cred_intervals
  
  sums <- c()
  for(i in 1:q){
    inti <- ocred[,,i]
    truei <- truemu[,i]
    indicatori <- (inti[,1] <= truei) & (inti[,2] >= truei)
    sums <- append(sums, sum(indicatori))
  }
  ocovmu <- 1/(N_subj*N_subj*q)*sum(sums)
  return(c(ocovmu))
}

getvcoverage <- function(N_samp, q, qobs, mods){
  data <- mods$data
  vstar <- data$Vstar[,(qobs + 1):q]
  
  opostsum <- extract_posterior_summaries(mods$omod, mods$osamps, 1:q)
  
  opred <- opostsum$pred_intervals
  
  sums <- c()
  for(i in (qobs+1):q){
    inti <- opred[,,i]
    truei <- vstar[,(i-qobs)]
    indicatori <- (inti[,1] <= truei) & (inti[,2] >= truei)
    sums <- append(sums, sum(indicatori))
  }
  ocovv <- 1/(N_subj*N_samp*(q-qobs))*sum(sums)
  return(c(ocovv))
}

# Default: N_subj=5, N_samp=10, q=10, q_obs=7
N_subj <- 5; N_samp <- 10; q <- 10; qobs <- 7;

rmses <- matrix(nrow = 1, ncol = length(simnames))
for (i in 1:length(simnames)){
  simname <- simnames[i]
  mods <- getsim_misspec_ours(simname)
  rmses[, i] <- getrmses(N_samp, q, qobs, mods)
}

rmspes <- matrix(nrow = 1, ncol = length(simnames))
for (i in 1:length(simnames)){
  simname <- simnames[i]
  mods <- getsim_misspec_ours(simname)
  rmspes[, i] <- getrmspes(N_samp, q, qobs, mods)
}

coveragemu <- matrix(nrow = 1, ncol = length(simnames))
for (i in 1:length(simnames)){
  simname <- simnames[i]
  mods <- getsim_misspec_ours(simname)
  coveragemu[, i] <- getmucoverage(N_samp, q, qobs, mods)
}

coveragev <- matrix(nrow = 1, ncol = length(simnames))
for (i in 1:length(simnames)){
  simname <- simnames[i]
  mods <- getsim_misspec_ours(simname)
  coveragev[, i] <- getvcoverage(N_samp, q, qobs, mods)
}

plot_rmses <- data.frame(rmse = t(rmses), simnames = simnames_short, 
                         Model = c(rep("Misspecified", times = 11), "Correctly specified"))

RMSEplot <- 
  plot_rmses %>% ggplot(aes(x = simnames, y = log(rmse), color = Model)) + 
  geom_point(cex = 4, alpha = 0.7, shape = 17) +
  scale_color_manual(values = c("#FC4E07", "blue")) +
  xlab("Data generation") + ylab("log(RMSE)") + ggtitle("RMSE") + 
  theme(legend.position = c(0.22, 0.86),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 9))

plot_rmspes <- data.frame(rmspe = t(rmspes), simnames = simnames_short, 
                          Model = c(rep("Misspecified", times = 11), "Correctly specified"))

RMSPEplot <- 
  plot_rmspes %>% ggplot(aes(x = simnames, y = rmspe, color = Model)) + 
  geom_point(cex = 4, alpha = 0.7, shape = 17) +
  scale_color_manual(values = c("#FC4E07", "blue")) +
  xlab("Data generation") + ylab("RMSPE") + ggtitle("RMSPE") + 
  theme(legend.position = 'none',
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 9))

plot_covmu <- data.frame(covmu = t(coveragemu), simnames = simnames_short, 
                         Model = c(rep("Misspecified", times = 11), "Correctly specified"))

covmuplot <- 
  plot_covmu %>% ggplot(aes(x = simnames, y = covmu, color = Model)) + 
  geom_point(cex = 4, alpha = 0.7, shape = 17) +
  scale_color_manual(values = c("#FC4E07", "blue")) +
  xlab("Data generation") + ylab("Coverage") + 
  ggtitle(expression(paste("Interval coverage for ", mu))) + 
  theme(legend.position = 'none',
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 9))

plot_covv <- data.frame(covv = t(coveragev), simnames = simnames_short, 
                        Model = c(rep("Misspecified", times = 11), "Correctly specified"))

covvplot <- 
  plot_covv %>% ggplot(aes(x = simnames, y = covv, color = Model)) + 
  geom_point(cex = 4, alpha = 0.7, shape = 17) +
  scale_color_manual(values = c("#FC4E07", "blue")) +
  xlab("Data generation") + ylab("Coverage") + 
  ggtitle("Interval coverage for V") + 
  theme(legend.position = 'none',
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 9))

ggsave("./sim-vis/figures/RMSE_misspec_all.png", RMSEplot, width = 6, height = 5)
ggsave("./sim-vis/figures/RMSPE_misspec_all.png", RMSPEplot, width = 6, height = 5)
ggsave("./sim-vis/figures/covmu_misspec_all.png", covmuplot, width = 6, height = 5)
ggsave("./sim-vis/figures/covv_misspec_all.png", covvplot, width = 6, height = 5)
