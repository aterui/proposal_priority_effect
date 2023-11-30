
# setup -------------------------------------------------------------------

source("code/library.R")
source("code/function.R")

df_sim <- readRDS("output/simulation.rds") %>% 
  mutate(stability = ifelse(eigen_max > 1, "unstable", "stable"))

df_plot <- df_sim %>% 
  drop_na(eigen_max)

# plot --------------------------------------------------------------------

g1 <- df_plot %>% 
  filter(sd_r == 0) %>% 
  ggplot(aes(x = stability,
             y = p)) +
  geom_violin(aes(fill = stability),
              alpha = 0.1,
              draw_quantiles = 0.5) +
  geom_jitter(height = 0,
              width = 0.1,
              alpha = 0.25,
              aes(color = eigen_max)) +
  facet_grid(rows = vars(n_timestep, factor_sd_a1),
             cols = vars(n_species),
             labeller = label_both)
  

g2 <- df_plot %>% 
  filter(sd_r == 0.1) %>% 
  ggplot(aes(x = stability,
             y = p)) +
  geom_violin(aes(fill = stability),
              alpha = 0.1,
              draw_quantiles = 0.5) +
  geom_jitter(height = 0,
              width = 0.1,
              alpha = 0.25) +
  facet_grid(rows = vars(n_timestep, factor_sd_a1),
             cols = vars(n_species),
             labeller = label_both)

g1 | g2


w1 <- df_plot %>% 
  filter(sd_r == 0) %>% 
  ggplot(aes(x = eigen_max,
             y = p,
             color = stability)) +
  geom_point(alpha = 0.25) +
  scale_x_continuous(trans = "log10") +
  facet_wrap(facets = ~ n_species + n_timestep + factor_sd_a1,
             #rows = vars(n_timestep, factor_sd_a1),
             #cols = vars(n_species),
             labeller = label_both,
             ncol = 2,
             nrow = 4,
             scales = "free_x")

# g1 <- df_plot %>%
#   filter(sd_r == 0) %>%
#   ggplot(aes(x = factor(factor_a1),
#              y = p,
#              color = stability)) +
#   geom_boxplot() +
#   geom_jitter(width = 0.1) +
#   facet_grid(rows = vars(n_timestep, factor_sd_a1),
#              cols = vars(n_species),
#              labeller = label_both) +
#   labs(x = "Relative strength of interspecific competition",
#        y = "Pr(b > b_0)")
# 
# g2 <- df_plot %>% 
#   filter(sd_r == 0.25) %>% 
#   ggplot(aes(x = factor(factor_a1),
#              y = p)) +
#   geom_boxplot() +
#   geom_jitter(width = 0.1) +
#   facet_grid(rows = vars(n_timestep, factor_sd_a1),
#              cols = vars(n_species),
#              labeller = label_both) +
#   labs(x = "Relative strength of interspecific competition",
#        y = "Pr(b > b_0)",
#        color = "SD in r")
