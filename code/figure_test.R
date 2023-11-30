
# setup -------------------------------------------------------------------

source("code/library.R")
source("code/function.R")

df_sim <- readRDS("output/simulation.rds") %>% 
  mutate(stability = ifelse(eigen_max > 1, "Sensitive", "Insensitive"))

df_plot <- df_sim %>% 
  drop_na(eigen_max)

# plot --------------------------------------------------------------------

label <- c('10' = 'Time series = 10',
           '30' = 'Time series = 30',
           '5' = "Species number = 5",
           '15' = "Species number = 15")

w <- df_plot %>% 
  ggplot(aes(x = eigen_max,
             y = p,
             color = stability)) +
  geom_point(alpha = 0.25,
             size = 0.5) +
  scale_x_continuous(trans = "log10") +
  facet_grid(rows = vars(n_timestep),
             cols = vars(n_species),
             labeller = as_labeller(label),
             scales = "free_x") +
  coord_cartesian(xlim = c(0.5, 10)) +
  labs(x = expression("Leading eigen value"~"|"*lambda[max]*"|"),
       y = expression("Pr("*delta[obs]~"<"~delta[null]*")"),
       color = "Sensitivity") +
  MetBrewer::scale_color_met_d("Hiroshige",
                               direction = -1) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom")

ggsave(w,
       filename = "output/figure_eigen.pdf",
       height = 3.5,
       width = 3.5)
