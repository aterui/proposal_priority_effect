
# setup -------------------------------------------------------------------
rm(list = ls())
source(here::here("code/library.R"))

#df_z <- readRDS(here::here("output/data_exponent_ssm.rds"))
df_z <- readRDS(here::here("output/data_exponent.rds"))

df_null <- df_z %>% 
  filter(neutral == 1) %>% 
  drop_na(z_dev) %>% 
  group_by(group) %>% 
  mutate(high = quantile(z_dev, 0.75),
         med = quantile(z_dev, 0.5),
         low = quantile(z_dev, 0.25)) %>% 
  ungroup()

# plot --------------------------------------------------------------------

ts <- 20
s <- 20
K <- 500

df_z_plot <- df_z %>%
  filter(nt == ts,
         nsp == s,
         min_k == K)

df_null_plot <- df_null %>%
  filter(nt == ts,
         nsp == s,
         min_k == K,
         sigma_alpha != 0)

df_z_plot %>% 
  filter(neutral == 0) %>% 
  ggplot(aes(x = factor(alpha),
             y = z_dev,
             color = factor(sigma_alpha),
             fill = factor(sigma_alpha))) +
  geom_jitter(alpha = 0.25,
              size = 0.5) +
  geom_boxplot(alpha = 0.1) +
  geom_hline(yintercept = quantile(df_null$z_dev,
                                   c(0.25, 0.5, 0.75)),
             linetype = "dashed",
             color = grey(0.3)) +
  scale_y_continuous(trans = "log10") +
  facet_grid(rows = vars(sigma_alpha),
             cols = vars(min_r, max_r),
             scales = "free") +
  labs(color = expression(sigma[alpha]),
       fill = expression(sigma[alpha])) +
  theme_bw() +
  theme(panel.grid = element_blank())

df_z_plot %>%
  filter(max_r != 2.5,
         sigma_alpha != 0.0001) %>%
  ggplot(aes(x = z_dev,
             color = factor(neutral),
             fill = factor(neutral))) +
  geom_density(alpha = 0.1) +
  scale_x_continuous(trans = "log10")
