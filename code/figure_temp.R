
# setup -------------------------------------------------------------------
rm(list = ls())
source(here::here("code/library.R"))

df_alpha <- readRDS(here::here("output/data_sim_ssm.rds")) %>% 
  group_by(group, replicate) %>% 
  summarize(alpha = unique(alpha))

df_z <- readRDS(here::here("output/data_exponent_ssm.rds")) %>% 
  left_join(df_alpha, by = c("group", "replicate"))

# plot --------------------------------------------------------------------

df_z %>%
  filter(nt == 20,
         nsp == 20,
         min_k == 1000,
         max_k == min_k,
         sigma_alpha == 0) %>% 
  ggplot(aes(x = alpha,
             y = z_dev,
             color = factor(sigma_alpha))) +
  geom_point(alpha = 0.1) +
  facet_grid(rows = vars(max_r),
             cols = vars(min_r),
             scales = "free") +
  geom_smooth(method = "lm") +
  #geom_hline(yintercept = 1) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(panel.grid = element_blank())
  
