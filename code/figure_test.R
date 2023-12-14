
# setup -------------------------------------------------------------------

source("code/library.R")
source("code/function.R")

df_sim <- readRDS("output/simulation.rds") %>% 
  mutate(stability = ifelse(eigen_max > 1, "Sensitive", "Insensitive"))

df_plot <- df_sim %>% 
  drop_na(eigen_max, p)

# df_neut <- df_sim %>% 
#   filter(sd_r == 0,
#          factor_sd_a1 == 0,
#          factor_a1 == 1) %>% 
#   group_by(n_timestep, n_species) %>% 
#   reframe(med = median(p),
#           up = quantile(p, 0.75),
#           low = quantile(p, 0.25))

# plot --------------------------------------------------------------------

label <- c('10' = 'Time series = 10',
           '30' = 'Time series = 30',
           '5' = "Species number = 5",
           '15' = "Species number = 15")

g_fill <- df_plot %>%
  ggplot(aes(fill = stability)) +
  geom_density(aes(x = p),
               color = "grey",
               alpha = 0.3,
               position = "fill") +
  facet_grid(rows = vars(n_timestep),
             cols = vars(n_species),
             labeller = as_labeller(label),
             scales = "free_x") +
  labs(x = expression("Pr("*delta[obs]~">"~delta[null]*")"),
       y = "Relative Frequency",
       fill = "Sensitivity") +
  MetBrewer::scale_color_met_d("Hiroshige",
                               direction = -1) +
  MetBrewer::scale_fill_met_d("Hiroshige",
                              direction = -1) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom")

ggsave(g_fill,
       filename = "output/figure_prop.pdf",
       height = 4.5,
       width = 4.5)

# boxplot option
g_box <- df_plot %>%
  # filter(a0 == 0.01,
  #        sd_r == 0) %>% 
  ggplot(aes(x = stability,
             y = p,
             color = stability,
             fill = stability)) +
  geom_jitter(height = 0,
              alpha = 0.3,
              size = 0.5) +
  geom_boxplot(alpha = 0.3,
               outlier.color = NA) +
  facet_grid(rows = vars(n_timestep),
             cols = vars(n_species),
             labeller = as_labeller(label),
             scales = "free_x") +
  labs(y = expression("Pr("*delta[obs]~">"~delta[null]*")"),
       x = "True sensitivity",
       color = "Sensitivity",
       fill = "Sensitivity") +
  guides(color = "none",
         fill = "none") +
  MetBrewer::scale_color_met_d("Hiroshige",
                               direction = -1) +
  MetBrewer::scale_fill_met_d("Hiroshige",
                              direction = -1) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom") +
  geom_ysidedensity(aes(x = stat(density)),
                    alpha = 0.50)

ggsave(g_box,
       filename = "output/figure_eigen.pdf",
       height = 4.5,
       width = 6)

# scatter plot option
# w <- df_plot %>%
#   ggplot(aes(x = eigen_max,
#              y = p,
#              color = stability)) +
#   geom_point(alpha = 0.25,
#              size = 1) +
#   scale_x_continuous(trans = "sqrt") +
#   facet_grid(rows = vars(n_timestep, factor_sd_a1),
#              cols = vars(n_species, sd_r),
#              labeller = label_both,#as_labeller(label),
#              scales = "free_x") +
#   coord_cartesian(xlim = c(0.5, 10)) +
#   labs(x = expression("Leading eigen value"~"|"*lambda[max]*"|"),
#        y = expression("Pr("*delta[obs]~"<"~delta[null]*")"),
#        color = "Sensitivity") +
#   MetBrewer::scale_color_met_d("Hiroshige",
#                                direction = -1) +
#   theme_bw() +
#   theme(strip.background = element_blank(),
#         panel.grid = element_blank(),
#         legend.position = "bottom")
# 
# ggsave(w,
#        filename = "output/figure_eigen.pdf",
#        height = 3.5,
#        width = 3.5)



# hsu results -------------------------------------------------------------

df_delta <- readRDS("output/simulation_hsu.rds")
df_hsu <- readRDS("data_fmt/data_hsu_fmt.rds")

df_x <- df_hsu %>% 
  group_by(replicate, treatment, light) %>% 
  summarize(x_max = max(x),
            y_max = max(y)) %>% 
  pivot_wider(id_cols = c("replicate", "light"),
              names_from = "treatment",
              values_from = c("x_max", "y_max")) %>% 
  ungroup() %>% 
  mutate(rr_x = abs(log(x_max_cp) - log(x_max_pc)),
         rr_y = abs(log(y_max_cp) - log(y_max_pc))) %>% 
  rowwise() %>% 
  summarize(replicate,
            light,
            log_rr = mean(rr_x, rr_y))

df_delta <- df_delta %>% 
  left_join(df_x, by = c("replicate", "light"))

g_hsu_box <- df_delta %>% 
  drop_na(p) %>% 
  mutate(sensitivity = ifelse(light == 0,
                              "Insensitive",
                              "Sensitive")) %>% 
  ggplot(aes(x = sensitivity,
             y = p,
             fill = sensitivity)) +
  geom_boxplot(outlier.color = NA,
              alpha = 0.5) +
  geom_jitter(height = 0,
             width = 0.1,
             alpha = 0.5) +
  MetBrewer::scale_fill_met_d("Hiroshige", direction = -1) +
  labs(x = "Sensitivity",
       y = expression("Pr("*delta[obs]~">"~delta[null]*")")) +
  guides(fill = "none") + 
  theme_bw() +
  theme(panel.grid = element_blank())

g_hsu <- df_delta %>% 
  drop_na(p) %>% 
  mutate(sensitivity = ifelse(light == 0,
                              "Insensitive",
                              "Sensitive")) %>% 
  ggplot(aes(x = p,
             y = log_rr,
             color = light)) +
  geom_point(size = 3) +
  MetBrewer::scale_color_met_c("Hiroshige", direction = -1) +
  labs(x = expression("Pr("*delta[obs]~">"~delta[null]*")"),
       y = "Strength of priority effects",
       color = "Light") +
  theme_bw() +
  theme(panel.grid = element_blank())


ggsave(g_hsu_box,
       filename = "output/figure_box_hsu.pdf",
       width = 4.5, height = 4)

ggsave(g_hsu,
       filename = "output/figure_scatter_hsu.pdf",
       width = 6, height = 4)