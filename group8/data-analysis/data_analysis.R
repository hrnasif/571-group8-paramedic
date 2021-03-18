library(paramedic)
library(rstan)
library(ggplot2)
library(dplyr)


data("group8_qpcr_data_v2")
data("group8_16S_data_v2")

group8mod <- paramedic::run_paramedic(W = group8_16S_data_v2, V = group8_qpcr_data_v2, n_samp = 10,
                         k = 0, n_iter = 16000, n_burnin = 14000, n_chains = 1)

saveRDS(group8mod, file = "group8/data-analysis/group8mod.rds")

new_W <- group8_16S_data_v2
new_V <- group8_qpcr_data_v2

new_W$subject_id <- 1:nrow(new_W)
new_V$subject_id <- 1:nrow(new_V)


# Previous model with shorter iterations
# stdmod <- paramedic::run_paramedic(W = new_W, V = new_V, n_samp = 1,
#                                    k = 0, n_iter = 16000, n_burnin = 14000, n_chains = 1)
# 
# saveRDS(stdmod, file = "group8/data-analysis/stdmod.rds")

stdmod <- paramedic::run_paramedic(W = new_W, V = new_V, n_samp = 1,
                                   k = 0, n_iter = 25000, n_burnin = 24000, n_chains = 1)

saveRDS(stdmod, file = "group8/data-analysis/stdmod_long.rds")


#### Obtain data for heatmaps
group8.samps <- rstan::extract(group8mod$stan_fit)
group8.post_summaries <- extract_posterior_summaries(group8mod$summary, group8.samps, 
                                                     taxa_of_interest = 1:27)

stdmod.samps <- rstan::extract(stdmod$stan_fit)
stdmod.post_summaries <- extract_posterior_summaries(stdmod$summary, stdmod.samps,
                                                     taxa_of_interest = 1:27)

grp8_min <- min(log(group8.post_summaries$estimates))
grp8_range <- max(log(group8.post_summaries$estimates)) - grp8_min
normalize_grp8 <- function(x){ (x - grp8_min)/(grp8_range) }

group8_log_estimates <- log(group8.post_summaries$estimates) %>% as.data.frame() %>%
  mutate_all(normalize_grp8) %>% cbind(group8_16S_data_v2[,1:2])



std_min <- min(log(stdmod.post_summaries$estimates))
std_range <- max(log(stdmod.post_summaries$estimates)) - std_min
normalize_std <- function(x){ (x - std_min)/(std_range) }

stdmod_log_estimates <- log(stdmod.post_summaries$estimates) %>% as.data.frame() %>% 
  mutate_all(normalize_std) %>% cbind(group8_16S_data_v2[,1:2])


colnames(group8_log_estimates)[1:27] <- colnames(group8_16S_data_v2)[3:29]
colnames(stdmod_log_estimates)[1:27] <- colnames(group8_16S_data_v2)[3:29]

group8_log_estimates_long <- tidyr::pivot_longer(group8_log_estimates,
                                                 cols = !ends_with("id"),
                                                 names_to = "taxon",
                                                 values_to = "log_mu")

stdmod_log_estimates_long <- tidyr::pivot_longer(stdmod_log_estimates,
                                                 cols = !ends_with("id"),
                                                 names_to = "taxon",
                                                 values_to = "log_mu")


### Make Heatmaps
group8.heatmap <- ggplot(data = group8_log_estimates_long, 
                         mapping = aes(x = interaction(as.factor(sample_id), as.factor(subject_id), lex.order = T), 
                                       y = reorder(taxon, log_mu), fill = log_mu)) +
  geom_tile() + 
  annotate(geom = "text", x = seq_len(nrow(group8_log_estimates)), y = -1, label = group8_log_estimates$sample_id, size = 3) +
  annotate(geom = "text", x = 5.5 + 10 * (0:3), y = -3, label = unique(group8_log_estimates$subject_id), size = 4) +
  coord_cartesian(ylim = c(-4, 28), expand = FALSE, clip = "off") +
  theme(
        axis.title.x = element_text(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "Estimated Concentration: with correlation") +
  xlab(label = "Subject and Sample IDs") + ylab(label = "Taxon") +
  guides(fill = guide_legend(title = latex2exp::TeX("Normalized $log(\\hat{\\mu}_{ijt})")))

group8.heatmap

stdmod.heatmap <- ggplot(data = stdmod_log_estimates_long, 
                         mapping = aes(x = interaction(as.factor(sample_id), as.factor(subject_id), lex.order = T), 
                                       y = reorder(taxon, log_mu), fill = log_mu)) +
  geom_tile() + 
  annotate(geom = "text", x = seq_len(nrow(group8_log_estimates)), y = -1, label = group8_log_estimates$sample_id, size = 3) +
  annotate(geom = "text", x = 5.5 + 10 * (0:3), y = -3, label = unique(group8_log_estimates$subject_id), size = 4) +
  coord_cartesian(ylim = c(-4, 28), expand = FALSE, clip = "off") +
  theme(
    axis.title.x = element_text(),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "Estimated Concentration: without correlation") +
  xlab(label = "Subject and Sample IDs") + ylab(label = "Taxon") +
  guides(fill = guide_legend(title = latex2exp::TeX("Normalized $log(\\hat{\\mu}_{ijt})")))

stdmod.heatmap


### Line Plot - by sample

qpcr_props <- prop.table(as.matrix(group8_qpcr_data_v2[-c(1,2)]), margin = 1)
qpcr_props <- cbind(as.matrix(group8_qpcr_data_v2[,c(1,2)]), qpcr_props) %>% as.data.frame() %>% 
  tidyr::pivot_longer(cols = !ends_with("id"), 
                      names_to = "taxon",
                      values_to = "qpcr_props")


group8mod_props <- prop.table(group8.post_summaries$estimates[,1:7], margin = 1) %>% 
  cbind(as.matrix(group8_qpcr_data_v2[,c(1,2)])) %>% as.data.frame() %>%
  setNames(colnames(group8_qpcr_data_v2)[c(3:9, 1,2)]) %>% 
  tidyr::pivot_longer(cols = !ends_with("id"),
                      names_to = "taxon",
                      values_to = "est_row_props") %>% mutate(Model = "group8")

stdmod_props <- prop.table(stdmod.post_summaries$estimates[,1:7], margin = 1) %>% 
  cbind(as.matrix(group8_qpcr_data_v2[,c(1,2)])) %>% as.data.frame() %>%
  setNames(colnames(group8_qpcr_data_v2)[c(3:9, 1,2)]) %>% 
  tidyr::pivot_longer(cols = !ends_with("id"),
                      names_to = "taxon",
                      values_to = "est_row_props") %>% mutate(Model = "Williamson")

model_props <- rbind(group8mod_props, stdmod_props)

full_props <- left_join(qpcr_props, model_props, by = c("subject_id", "sample_id", "taxon"))
full_props$Model <- as.factor(full_props$Model)

props_line_sample <- ggplot(data = full_props, aes(x = qpcr_props, y = est_row_props)) + 
  geom_point(aes(shape = Model, color = Model), size = 2) + 
  geom_abline(intercept = 0, slope = 1, color = "black") +
  labs(title = "Estimated vs. Observed Relative Abundance by Sample", x = "Observed (qPCR)",
       y = "Estimated") + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
props_line_sample
ggsave("group8/data-analysis/sample_lineplot_long.png", plot = props_line_sample, width = 8, height = 3.5)



### Line plot - by subject

qpcr_props_subj <- group8_qpcr_data_v2 %>% group_by(subject_id) %>% select(!sample_id) %>% summarise_all(sum)

qpcr_props_subj <- prop.table(as.matrix(qpcr_props_subj[,-1]), margin = 1)
qpcr_props_subj <- cbind(as.matrix(unique(group8_qpcr_data_v2[,1])), qpcr_props_subj) %>% as.data.frame() %>% 
  tidyr::pivot_longer(cols = !ends_with("id"), 
                      names_to = "taxon",
                      values_to = "qpcr_props")


group8mod_props_subj <- cbind(as.matrix(group8_qpcr_data_v2[,1]), group8.post_summaries$estimates[,1:7]) %>% 
  as.data.frame %>% group_by(subject_id) %>% summarise_all(sum)

group8mod_props_subj <- prop.table(as.matrix(group8mod_props_subj[,-1]), margin = 1) %>% 
  cbind(as.matrix(as.matrix(unique(group8_qpcr_data_v2[,1])))) %>% as.data.frame() %>%
  setNames(colnames(group8_qpcr_data_v2)[c(3:9,1)]) %>% 
  tidyr::pivot_longer(cols = !ends_with("id"),
                      names_to = "taxon",
                      values_to = "est_row_props") %>% mutate(Model = "group8")



stdmod_props_subj <- cbind(as.matrix(group8_qpcr_data_v2[,1]), stdmod.post_summaries$estimates[,1:7]) %>% 
  as.data.frame %>% group_by(subject_id) %>% summarise_all(sum)

stdmod_props_subj <- prop.table(as.matrix(stdmod_props_subj[,-1]), margin = 1) %>% 
  cbind(as.matrix(as.matrix(unique(group8_qpcr_data_v2[,1])))) %>% as.data.frame() %>%
  setNames(colnames(group8_qpcr_data_v2)[c(3:9,1)]) %>% 
  tidyr::pivot_longer(cols = !ends_with("id"),
                      names_to = "taxon",
                      values_to = "est_row_props") %>% mutate(Model = "Williamson")

model_props_subj <- rbind(group8mod_props_subj, stdmod_props_subj)

full_props_subj <- left_join(qpcr_props_subj, model_props_subj, by = c("subject_id", "taxon"))
full_props_subj$Model <- as.factor(full_props_subj$Model)
full_props_subj$subject_id <- as.factor(full_props_subj$subject_id)

props_line_sample_subj <- ggplot(data = full_props_subj, aes(x = qpcr_props, y = est_row_props)) + 
  geom_point(aes(shape = Model, color = subject_id), size = 3) + 
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "Estimated vs. Observed Relative Abundance by Subject", x = "Observed (qPCR)",
       y = "Estimated") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


props_line_sample_subj
ggsave("group8/data-analysis/subj_lineplot_long.png", plot = props_line_sample_subj, width = 8, height = 3.5)


