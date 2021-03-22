library(ggplot2)
setwd("./group8/")
source("./sim-vis/extractinits.R")
source("./sim-vis/simfns.R")

# Vary N_samp, q=10, q_obs=7, (N_subj=5)
N_subj <- 5
q <- 10
qobs <- 7
N_samp <- c(5, 10, 15, 20)


rmses <- matrix(nrow = 3, ncol = length(N_samp))
for (i in 1:length(N_samp)){
  nsamp <- N_samp[i]
  rmses[, i] <- getrmses(nsamp, q, qobs)
}

rmspes <- matrix(nrow = 3, ncol = length(N_samp))
for (i in 1:length(N_samp)){
  nsamp <- N_samp[i]
  rmspes[, i] <- getrmspes(nsamp, q, qobs)
}

coveragemu <- matrix(nrow = 2, ncol = length(N_samp))
for (i in 1:length(N_samp)){
  nsamp <- N_samp[i]
  coveragemu[, i] <- getmucoverage(nsamp, q, qobs)
}

coveragev <- matrix(nrow = 2, ncol = length(N_samp))
for (i in 1:length(N_samp)){
  nsamp <- N_samp[i]
  coveragev[, i] <- getvcoverage(nsamp, q, qobs)
}

plot_rmses <- data.frame(rmse = c(t(rmses)),
                          N_samp = rep(N_samp, 3), 
                          Estimator = factor(rep(c("Proposed", "Williamson et al.", "Naive"), each = 4), 
                                             levels = c("Naive", "Williamson et al.", "Proposed")
                                             )
                         )
RMSEplot <- plot_rmses %>% ggplot(aes(x = N_samp, y = rmse, color = Estimator, shape = Estimator)) + 
  geom_point(cex = 4, alpha = 0.7) +
  scale_color_manual(values = c("#E7B800", "#2E9FDF", "#FC4E07")) +
  xlab("Number of samples per individual (m)") + ylab("RMSE") + ggtitle("RMSE") +
  theme(legend.position = c(0.22, 0.86),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15)
        #panel.grid.minor = element_line(color = "black")
        )

plot_rmspes <- data.frame(rmspe = c(t(rmspes)),
                         N_samp = rep(N_samp, 3), 
                         Estimator = factor(rep(c("Proposed", "Williamson et al.", "Naive"), each = 4), 
                                            levels = c("Naive", "Williamson et al.", "Proposed")
                                            )
                         )
RMSPEplot <- plot_rmspes %>% ggplot(aes(x = N_samp, y = rmspe, color = Estimator, shape = Estimator)) + 
              geom_point(cex = 4, alpha = 0.7) +
              scale_color_manual(values = c("#E7B800", "#2E9FDF", "#FC4E07")) +
              xlab("Number of samples per individual (m)") + ylab("RMSPE") + ggtitle("RMSPE") + 
              theme(legend.position = 'none',
                    axis.line = element_line(color = "black"),
                    panel.background = element_rect(fill = "white"),
                    text = element_text(size = 15))


plot_covmu <- data.frame(cov = c(t(coveragemu)),
                          N_samp = rep(N_samp, 2), 
                          Estimator = factor(rep(c("Proposed", "Williamson et al."), each = 4), 
                                             levels = c("Williamson et al.", "Proposed")
                                             )
                          )
covmuplot <- plot_covmu %>% ggplot(aes(x = N_samp, y = cov, color = Estimator, shape = Estimator)) + 
  geom_point(cex = 4, alpha = 0.7) +
  scale_color_manual(values = c("#2E9FDF", "#FC4E07")) +
  xlab("Number of samples per individual (m)") + ylab("Coverage") + 
  ggtitle(expression(paste("Interval coverage for ", mu))) + 
  theme(legend.position = 'none',
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15))

plot_covv <- data.frame(cov = c(t(coveragev)),
                         N_samp = rep(N_samp, 2), 
                         Estimator = factor(rep(c("Proposed", "Williamson et al."), each = 4), 
                                            levels = c("Williamson et al.", "Proposed")
                                            )
                         )
covvplot <- plot_covv %>% ggplot(aes(x = N_samp, y = cov, color = Estimator, shape = Estimator)) + 
  geom_point(cex = 4, alpha = 0.7) +
  scale_color_manual(values = c("#2E9FDF", "#FC4E07")) +
  xlab("Number of samples per individual (m)") + ylab("Coverage") + ggtitle("Interval coverage for V")+ 
  theme(legend.position = 'none',
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15))

ggsave("./sim-vis/figures/RMSE_Nsamp.png", RMSEplot, width = 5, height = 5)
ggsave("./sim-vis/figures/RMSPE_Nsamp.png", RMSPEplot, width = 5, height = 5)
ggsave("./sim-vis/figures/covmu_Nsamp.png", covmuplot, width = 5, height = 5)
ggsave("./sim-vis/figures/covv_Nsamp.png", covvplot, width = 5, height = 5)
