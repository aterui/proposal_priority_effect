
# setup -------------------------------------------------------------------
rm(list = ls())
source(here::here("code/library.R"))

df_alpha <- readRDS(here::here("output/data_sim_ssm.rds")) %>% 
  group_by(group, replicate) %>% 
  summarize(alpha = unique(alpha))

df_z <- readRDS(here::here("output/data_exponent_ssm.rds")) %>% 
  left_join(df_alpha, by = c("group", "replicate"))

df_null <- df_z %>% 
  filter(neutral == 1) %>% 
  group_by(group) %>% 
  mutate(up = quantile(z_dev, 0.75),
         z_dev = quantile(z_dev, 0.5),
         low = quantile(z_dev, 0.25),
         alpha = seq(0, 1, length = max(replicate)))

# plot --------------------------------------------------------------------

ts <- 20
s <- 20
K <- 1000

df_z_plot <- df_z %>%
  filter(nt == ts,
         nsp == s,
         min_k == K)

df_null_plot <- df_null %>%
  filter(nt == ts,
         nsp == s,
         min_k == K)

df_z_plot %>% 
  filter(neutral == 0,
         sigma_alpha == 0.25) %>% 
  ggplot(aes(x = alpha,
             y = z_dev,
             color = factor(sigma_alpha),
             fill = factor(sigma_alpha))) +
  geom_point(alpha = 0.2) +
  facet_grid(rows = vars(max_r),
             cols = vars(min_r),
             scales = "free") +
  geom_ribbon(data = df_null_plot,
              aes(ymin = low,
                  ymax = up),
              col = NA,
              fill = grey(0),
              alpha = 0.15) +
  geom_smooth(method = "lm") +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(panel.grid = element_blank())

# df_z_plot %>% 
#   filter(max_r != 2.5) %>% 
#   ggplot(aes(x = z_dev,
#              color = factor(neutral),
#              fill = factor(neutral))) +
#   geom_density(alpha = 0.1) +
#   scale_x_continuous(trans = "log10")
