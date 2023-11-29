
# setup -------------------------------------------------------------------

source("code/library.R")
source("code/function.R")

df_plot <- readRDS("output/simulation.rds")


# plot --------------------------------------------------------------------

df_plot %>% 
  ggplot(aes(x = factor(lambda_a1),
             y = p,
             color = factor(sd_r))) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  facet_grid(rows = vars(n_timestep, lambda_sd_a1),
             cols = vars(n_species),
             labeller = label_both) +
  labs(x = "Relative strength of interspecific competition",
       y = "Pr(b > b_0)",
       color = "SD in r")
