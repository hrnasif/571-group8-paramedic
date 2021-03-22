library(ggplot2)
setwd("./group8/")
source("./sim-vis/extractinits.R")
source("./sim-vis/simfns.R")

# Vary q_obs, q=10, N_samp=10 (defaults)
N_subj <- 5
q <- 10
qobs <- c(3, 5, 7)
N_samp <- 10


rmses <- matrix(nrow = 3, ncol = length(qobs))
for (i in 1:length(qobs)){
  q_obs <- qobs[i]
  mods <- getsim(N_samp, q, q_obs)
  rmses[, i] <- getrmses(N_samp, q, q_obs, mods)
}

rmspes <- matrix(nrow = 3, ncol = length(qobs))
for (i in 1:length(qobs)){
  q_obs <- qobs[i]
  mods <- getsim(N_samp, q, q_obs)
  rmspes[, i] <- getrmspes(N_samp, q, q_obs, mods)
}

coveragemu <- matrix(nrow = 2, ncol = length(qobs))
for (i in 1:length(qobs)){
  q_obs <- qobs[i]
  mods <- getsim(N_samp, q, q_obs)
  coveragemu[, i] <- getmucoverage(N_samp, q, q_obs, mods)
}

coveragev <- matrix(nrow = 2, ncol = length(qobs))
for (i in 1:length(qobs)){
  q_obs <- qobs[i]
  mods <- getsim(N_samp, q, q_obs)
  coveragev[, i] <- getvcoverage(N_samp, q, q_obs, mods)
}

plot_rmses <- data.frame(rmse = c(t(rmses)),
                         qobs = rep(qobs, 3), 
                         Estimator = factor(rep(c("Proposed", "Williamson et al.", "Naive"), each = length(qobs)), 
                                            levels = c("Naive", "Williamson et al.", "Proposed")
                         )
)
RMSEplot <- plot_rmses %>% ggplot(aes(x = qobs, y = rmse, color = Estimator, shape = Estimator)) + 
  geom_point(cex = 4, alpha = 0.7) +
  scale_color_manual(values = c("#E7B800", "#2E9FDF", "#FC4E07")) +
  xlab(expression('q'["obs"])) + ylab("RMSE") + ggtitle("RMSE") +
  theme(legend.position = c(0.22, 0.86),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15)
        #panel.grid.minor = element_line(color = "black")
  )

plot_rmspes <- data.frame(rmspe = c(t(rmspes)),
                          qobs = rep(qobs, 3), 
                          Estimator = factor(rep(c("Proposed", "Williamson et al.", "Naive"), each = length(qobs)), 
                                             levels = c("Naive", "Williamson et al.", "Proposed")
                          )
)
RMSPEplot <- plot_rmspes %>% ggplot(aes(x = qobs, y = rmspe, color = Estimator, shape = Estimator)) + 
  geom_point(cex = 4, alpha = 0.7) +
  scale_color_manual(values = c("#E7B800", "#2E9FDF", "#FC4E07")) +
  xlab(expression('q'["obs"])) + ylab("RMSPE") + ggtitle("RMSPE") + 
  theme(legend.position = 'none',
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15))


plot_covmu <- data.frame(cov = c(t(coveragemu)),
                         qobs = rep(qobs, 2), 
                         Estimator = factor(rep(c("Proposed", "Williamson et al."), each = length(qobs)), 
                                            levels = c("Williamson et al.", "Proposed")
                         )
)
covmuplot <- plot_covmu %>% ggplot(aes(x = qobs, y = cov, color = Estimator, shape = Estimator)) + 
  geom_point(cex = 4, alpha = 0.7) +
  scale_color_manual(values = c("#2E9FDF", "#FC4E07")) +
  xlab(expression('q'["obs"])) + ylab("Coverage") + 
  ggtitle(expression(paste("Interval coverage for ", mu))) + 
  theme(legend.position = 'none',
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15))

plot_covv <- data.frame(cov = c(t(coveragev)),
                        qobs = rep(qobs, 2), 
                        Estimator = factor(rep(c("Proposed", "Williamson et al."), each = length(qobs)), 
                                           levels = c("Williamson et al.", "Proposed")
                        )
)
covvplot <- plot_covv %>% ggplot(aes(x = qobs, y = cov, color = Estimator, shape = Estimator)) + 
  geom_point(cex = 4, alpha = 0.7) +
  scale_color_manual(values = c("#2E9FDF", "#FC4E07")) +
  xlab(expression('q'["obs"])) + ylab("Coverage") + ggtitle("Interval coverage for V")+ 
  theme(legend.position = 'none',
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15))

ggsave("./sim-vis/figures/RMSE_qobs.png", RMSEplot, width = 5, height = 5)
ggsave("./sim-vis/figures/RMSPE_qobs.png", RMSPEplot, width = 5, height = 5)
ggsave("./sim-vis/figures/covmu_qobs.png", covmuplot, width = 5, height = 5)
ggsave("./sim-vis/figures/covv_qobs.png", covvplot, width = 5, height = 5)