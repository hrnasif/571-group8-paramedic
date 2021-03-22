library(ggplot2)
setwd("./group8/")
source("./sim-vis/extractinits.R")
source("./sim-vis/simfns.R")

# Vary q_obs, q=10, N_samp=10 (defaults)
N_subj <- 5
q <- c(5, 10, 15)
qobs <- 7
N_samp <- 10


rmses <- matrix(nrow = 3, ncol = length(q))
for (i in 1:length(q)){
  q_t <- q[i]
  if(i == 1){
    mods <- getsim(N_samp, q_t, 3)
    rmses[, i] <- getrmses(N_samp, q_t, 3, mods)
  }else{
    mods <- getsim(N_samp, q_t, qobs)
    rmses[, i] <- getrmses(N_samp, q_t, qobs, mods)
  }
}

rmspes <- matrix(nrow = 3, ncol = length(q))
for (i in 1:length(q)){
  q_t <- q[i]
  if(i == 1){
    mods <- getsim(N_samp, q_t, 3)
    rmspes[, i] <- getrmspes(N_samp, q_t, 3, mods)
  } else{
    mods <- getsim(N_samp, q_t, qobs)
    rmspes[, i] <- getrmspes(N_samp, q_t, qobs, mods)
  }
}

coveragemu <- matrix(nrow = 2, ncol = length(q))
for (i in 1:length(q)){
  q_t <- q[i]
  if(i == 1){
    mods <- getsim(N_samp, q_t, 3)
    coveragemu[, i] <- getmucoverage(N_samp, q_t, 3, mods)
  } else{
    mods <- getsim(N_samp, q_t, qobs)
    coveragemu[, i] <- getmucoverage(N_samp, q_t, qobs, mods)
  }
}

coveragev <- matrix(nrow = 2, ncol = length(q))
for (i in 1:length(q)){
  q_t <- q[i]
  if(i == 1){
    mods <- getsim(N_samp, q_t, 3)
    coveragev[, i] <- getvcoverage(N_samp, q_t, 3, mods)
  } else{
    mods <- getsim(N_samp, q_t, qobs)
    coveragev[, i] <- getvcoverage(N_samp, q_t, qobs, mods)
  }
}

plot_rmses <- data.frame(rmse = c(t(rmses)),
                         q = rep(q, 3), 
                         Estimator = factor(rep(c("Proposed", "Williamson et al.", "Naive"), each = length(q)), 
                                            levels = c("Naive", "Williamson et al.", "Proposed")
                         )
)
RMSEplot <- plot_rmses %>% ggplot(aes(x = q, y = rmse, color = Estimator, shape = Estimator)) + 
  geom_point(cex = 4, alpha = 0.7) +
  scale_color_manual(values = c("#E7B800", "#2E9FDF", "#FC4E07")) +
  xlab("q") + ylab("RMSE") + ggtitle("RMSE") +
  theme(legend.position = c(0.22, 0.86),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15)
        #panel.grid.minor = element_line(color = "black")
  )

plot_rmspes <- data.frame(rmspe = c(t(rmspes)),
                          q = rep(q, 3), 
                          Estimator = factor(rep(c("Proposed", "Williamson et al.", "Naive"), each = length(q)), 
                                             levels = c("Naive", "Williamson et al.", "Proposed")
                          )
)
RMSPEplot <- plot_rmspes %>% ggplot(aes(x = q, y = rmspe, color = Estimator, shape = Estimator)) + 
  geom_point(cex = 4, alpha = 0.7) +
  scale_color_manual(values = c("#E7B800", "#2E9FDF", "#FC4E07")) +
  xlab("q") + ylab("RMSPE") + ggtitle("RMSPE") + 
  theme(legend.position = 'none',
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15))


plot_covmu <- data.frame(cov = c(t(coveragemu)),
                         q = rep(q, 2), 
                         Estimator = factor(rep(c("Proposed", "Williamson et al."), each = length(q)), 
                                            levels = c("Williamson et al.", "Proposed")
                         )
)
covmuplot <- plot_covmu %>% ggplot(aes(x = q, y = cov, color = Estimator, shape = Estimator)) + 
  geom_point(cex = 4, alpha = 0.7) +
  scale_color_manual(values = c("#2E9FDF", "#FC4E07")) +
  xlab("q") + ylab("Coverage") + 
  ggtitle(expression(paste("Interval coverage for ", mu))) + 
  theme(legend.position = 'none',
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15))

plot_covv <- data.frame(cov = c(t(coveragev)),
                        q = rep(q, 2), 
                        Estimator = factor(rep(c("Proposed", "Williamson et al."), each = length(q)), 
                                           levels = c("Williamson et al.", "Proposed")
                        )
)
covvplot <- plot_covv %>% ggplot(aes(x = q, y = cov, color = Estimator, shape = Estimator)) + 
  geom_point(cex = 4, alpha = 0.7) +
  scale_color_manual(values = c("#2E9FDF", "#FC4E07")) +
  xlab("q") + ylab("Coverage") + ggtitle("Interval coverage for V")+ 
  theme(legend.position = 'none',
        axis.line = element_line(color = "black"),
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 15))

ggsave("./sim-vis/figures/RMSE_q.png", RMSEplot, width = 5, height = 5)
ggsave("./sim-vis/figures/RMSPE_q.png", RMSPEplot, width = 5, height = 5)
ggsave("./sim-vis/figures/covmu_q.png", covmuplot, width = 5, height = 5)
ggsave("./sim-vis/figures/covv_q.png", covvplot, width = 5, height = 5)
