#setwd("./group8/")
source("./sim-vis/simfns.R")
source("./sim-vis/extractinits.R")
simnames <- c("misspec-gamma-gamma-poisson-mult", "misspec-normal-gamma-poisson-mult", "misspec-gamma-normal-poisson-mult")
simnames_short <- c("gamma, gamma", "normal, gamma", "gamma, normal")

# Vary sim-name (misspecify model)
# Default: N_subj=5, N_samp=10, q=15, q_obs=7
N_subj <- 5; N_samp <- 10; q <- 10; qobs <- 7;

rmses <- matrix(nrow = 3, ncol = length(simnames))
for (i in 1:length(simnames)){
  simname <- simnames[i]
  mods <- getsim_misspec(simname)
  rmses[, i] <- getrmses(N_samp, q, qobs, mods)
}

rmspes <- matrix(nrow = 3, ncol = length(simnames))
for (i in 1:length(simnames)){
  simname <- simnames[i]
  mods <- getsim_misspec(simname)
  rmspes[, i] <- getrmspes(N_samp, q, qobs, mods)
}

coveragemu <- matrix(nrow = 2, ncol = length(simnames))
for (i in 1:length(simnames)){
  simname <- simnames[i]
  mods <- getsim_misspec(simname)
  coveragemu[, i] <- getmucoverage(N_samp, q, qobs, mods)
}

coveragev <- matrix(nrow = 2, ncol = length(simnames))
for (i in 1:length(simnames)){
  simname <- simnames[i]
  mods <- getsim_misspec(simname)
  coveragev[, i] <- getvcoverage(N_samp, q, qobs, mods)
}

# plot_rmses <- data.frame(rmse = c(t(rmses)),
#                          simnames = factor(rep(simnames_short, 3)), 
#                          Estimator = factor(rep(c("Proposed", "Williamson et al.", "Naive"), each = length(simnames)), 
#                                             levels = c("Naive", "Williamson et al.", "Proposed")
#                          )
# )
plot_rmses <- data.frame(rmse = c(t(rmses[1:2,])),
                         simnames = factor(rep(simnames_short, 2)), 
                         Estimator = factor(rep(c("Proposed", "Williamson et al."), each = length(simnames)), 
                                            levels = c("Williamson et al.", "Proposed")
                         )
)
RMSEplot <- plot_rmses %>% ggplot(aes(x = simnames, y = rmse, color = Estimator, shape = Estimator)) + 
  geom_point(cex = 4, alpha = 0.7) +
  scale_color_manual(values = c("#2E9FDF", "#FC4E07")) +
  scale_shape_manual(values = c(17, 15)) +
  xlab(expression(paste(mu, ", e data generation"))) + ylab("RMSE") + ggtitle("RMSE") +
  theme(legend.position = c(0.22, 0.5),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15)
        #panel.grid.minor = element_line(color = "black")
  )

plot_rmspes <- data.frame(rmspe = c(t(rmspes[1:2,])),
                          simnames = rep(simnames_short, 2), 
                          Estimator = factor(rep(c("Proposed", "Williamson et al."), each = length(simnames)), 
                                             levels = c("Williamson et al.", "Proposed")
                          )
)
RMSPEplot <- plot_rmspes %>% ggplot(aes(x = simnames, y = rmspe, color = Estimator, shape = Estimator)) + 
  geom_point(cex = 4, alpha = 0.7) +
  scale_color_manual(values = c("#2E9FDF", "#FC4E07")) +
  scale_shape_manual(values = c(17, 15)) +
  xlab(expression(paste(mu, ", e data generation"))) + ylab("RMSPE") + ggtitle("RMSPE") + 
  theme(legend.position = 'none',
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15))


plot_covmu <- data.frame(cov = c(t(coveragemu)),
                         simnames = rep(simnames_short, 2), 
                         Estimator = factor(rep(c("Proposed", "Williamson et al."), each = length(simnames)), 
                                            levels = c("Williamson et al.", "Proposed")
                         )
)
covmuplot <- plot_covmu %>% ggplot(aes(x = simnames, y = cov, color = Estimator, shape = Estimator)) + 
  geom_point(cex = 4, alpha = 0.7) +
  scale_color_manual(values = c("#2E9FDF", "#FC4E07")) +
  scale_shape_manual(values = c(17, 15)) +
  xlab(expression(paste(mu, ", e data generation"))) + ylab("Coverage") + 
  ggtitle(expression(paste("Interval coverage for ", mu))) + 
  theme(legend.position = 'none',
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15))

plot_covv <- data.frame(cov = c(t(coveragev)),
                        simnames = rep(simnames_short, 2), 
                        Estimator = factor(rep(c("Proposed", "Williamson et al."), each = length(simnames)), 
                                           levels = c("Williamson et al.", "Proposed")
                        )
)
covvplot <- plot_covv %>% ggplot(aes(x = simnames, y = cov, color = Estimator, shape = Estimator)) + 
  geom_point(cex = 4, alpha = 0.7) +
  scale_color_manual(values = c("#2E9FDF", "#FC4E07")) +
  scale_shape_manual(values = c(17, 15)) +
  xlab(expression(paste(mu, ", e data generation"))) + ylab("Coverage") + ggtitle("Interval coverage for V")+ 
  theme(legend.position = 'none',
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15))

ggsave("./sim-vis/figures/RMSE_misspec.png", RMSEplot, width = 5, height = 5)
ggsave("./sim-vis/figures/RMSPE_misspec.png", RMSPEplot, width = 5, height = 5)
ggsave("./sim-vis/figures/covmu_misspec.png", covmuplot, width = 5, height = 5)
ggsave("./sim-vis/figures/covv_misspec.png", covvplot, width = 5, height = 5)
