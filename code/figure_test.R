
# setup -------------------------------------------------------------------

source("code/library.R")
source("code/function.R")

# simulation plot ---------------------------------------------------------

## read data
df_sim <- readRDS("output/simulation_stan.rds") %>% 
  mutate(stability = ifelse(eigen_max > 1, "Sensitive", "Insensitive"))

df_plot <- df_sim %>% 
  drop_na(eigen_max, p)

## label
label_t <- c('10' = 'Time series = 10',
             '30' = 'Time series = 30')

label_s <- c('2' = "Species number = 2",
             '5' = "Species number = 5",
             '10' = "Species number = 10")

## filled frequencies
g_fill <- df_plot %>%
  ggplot(aes(fill = stability)) +
  geom_density(aes(x = p),
               color = "grey",
               alpha = 0.3,
               position = "fill") +
  facet_grid(rows = vars(n_timestep),
             cols = vars(n_species),
             labeller = labeller(n_timestep = label_t,
                                 n_species = label_s),
             scales = "free_x") +
  labs(x = expression("Competitive exceedance"~Psi),
       y = "Relative frequency",
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
       width = 6)

## boxplot
g_box <- df_plot %>%
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
             labeller = labeller(n_timestep = label_t,
                                 n_species = label_s),
             scales = "free_x") +
  labs(y = expression("Competitive exceedance"~Psi),
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
       filename = "output/figure_eigen_box.pdf",
       height = 4.5,
       width = 6)

## scatter plot option
g_scatter <- df_plot %>%
  ggplot(aes(x = eigen_max,
             y = p,
             color = stability)) +
  geom_point(alpha = 0.4,
             size = 0.25) +
  #scale_x_continuous(trans = "sqrt") +
  #coord_cartesian(xlim = c(0, .5)) +
  facet_grid(rows = vars(n_timestep),
             cols = vars(n_species),
             labeller = labeller(n_timestep = label_t,
                                 n_species = label_s),
             scales = "free") +
  labs(x = expression("Leading eigen value"~"|"*lambda[max]*"|"),
       y = expression("Competitive exceedance"~Psi),
       color = "Sensitivity",
       fill = "Sensitivity") +
  geom_vline(xintercept = 1,
             linewidth = 0.25,
             linetype = "dashed",
             alpha = 0.4) +
  geom_hline(yintercept = 0.5,
             linewidth = 0.25,
             linetype = "dashed",
             alpha = 0.4) +
  geom_ysidedensity(aes(x = stat(density),
                        fill = stability),
                    alpha = 0.50,
                    position = "fill") +
  MetBrewer::scale_color_met_d("Hiroshige",
                               direction = -1) +
  MetBrewer::scale_fill_met_d("Hiroshige",
                              direction = -1) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        ggside.axis.text = element_blank(),
        ggside.axis.ticks = element_blank())
  
ggsave(g_scatter,
       filename = "output/figure_eigen_scatter.pdf",
       height = 5,
       width = 7)


# hsu results -------------------------------------------------------------

## read and format data
df_delta_hsu <- readRDS("output/simulation_hsu.rds") %>% 
  mutate(source = "Hsu and Moeller 2021")
df_delta_ojima <- readRDS("output/simulation_ojima.rds") %>% 
  mutate(source = "Ojima et al. 2022")
df_hsu <- readRDS("data_fmt/data_hsu_fmt.rds")

df_delta <- df_delta_hsu %>% 
  bind_rows(df_delta_ojima) %>% 
  select(source, delta, nsim, n, p, n_sample, p_neg, log_rr) %>% 
  drop_na(p)

## hsu - box plot
g_hsu_box <- df_delta_hsu %>% 
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
       y = expression("Competitive exceedance"~Psi)) +
  guides(fill = "none") + 
  theme_bw() +
  theme(panel.grid = element_blank())

## hsu - scatter plot 1
g_hsu_light <- df_delta_hsu %>% 
  drop_na(p) %>% 
  ggplot(aes(x = p,
             y = log_rr,
             color = light)) +
  geom_point(size = 3) +
  #MetBrewer::scale_color_met_d("Hiroshige", direction = -1) +
  MetBrewer::scale_color_met_c("Hiroshige", direction = -1) +
  labs(x = "Competitive exceedance Q",#expression("Pr("*delta[obs]~">"~delta[null]*")"),
       y = "Strength of priority effects",
       color = "Light") +
  theme_bw() +
  theme(panel.grid = element_blank())

## hsu - scatter plot 2
g_hsu_pr <- df_delta_hsu %>% 
  drop_na(p) %>% 
  ggplot(aes(x = p,
             y = log_rr)) +
  geom_point(size = 3) +
  #MetBrewer::scale_color_met_d("Hiroshige", direction = -1) +
  MetBrewer::scale_color_met_c("Hiroshige", direction = -1) +
  labs(x = expression("Competitive exceedance"~Psi),
       y = "Strength of priority effects",
       color = "Pr(negative)") +
  theme_bw() +
  theme(panel.grid = element_blank())

## combined - scatter plot
g_merge <- df_delta %>% 
  ggplot(aes(x = p,
             y = log_rr)) +
  geom_point(size = 1.5,
             alpha = 0.5,
             shape = "circle",
             data = . %>% filter(source == "Hsu and Moeller 2021")) +
  geom_point(size = 1.5,
             alpha = 0.5,
             shape = "triangle",
             data = . %>% filter(source == "Ojima et al. 2022")) +
  geom_smooth(color = grey(0.25),
              method = "lm",
              alpha = 0.3) +
  MetBrewer::scale_color_met_d("Hiroshige", direction = -1) +
  #MetBrewer::scale_color_met_c("Hiroshige", direction = -1) +
  labs(x = expression("Competitive exceedance"~Psi),
       y = "Strength of priority effects") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous(trans = "log10") +
  guides(shape = "none")

## export
# ggsave(g_hsu_box,
#        filename = "output/figure_hsu_box.pdf",
#        width = 4.5, height = 4)
# 
# ggsave(g_hsu_light,
#        filename = "output/figure_hsu_scatter_light.pdf",
#        width = 6, height = 4)
# 
# ggsave(g_hsu_pr,
#        filename = "output/figure_hsu_scatter_pr.pdf",
#        width = 6, height = 4)

ggsave(g_merge,
       filename = "output/figure_exp.pdf",
       width = 4, height = 3.5)

# help plot: experiment ---------------------------------------------------

# ## time-series
# df_t <- df_hsu %>% 
#   filter(n > 0) %>% 
#   filter(predictor == "x") %>%
#   select(-species) %>%
#   group_by(replicate, treatment, light) %>% 
#   mutate(t = day - min(day0) + 1) %>% 
#   pivot_longer(cols = c("x", "y"),
#                names_to = "species",
#                values_to = "density")
# 
# df_t %>% 
#   ggplot(aes(x = t,
#              y = density,
#              color = species)) +
#   geom_point() +
#   geom_line() +
#   scale_y_continuous(trans = "sqrt") +
#   facet_grid(rows = vars(replicate, treatment),
#              cols = vars(light))
# 
# 
# 
# ## community level plot
# df_n <- df_hsu %>% 
#   #filter(x > 0 | y > 0) %>% 
#   group_by(day, replicate, treatment, light) %>% 
#   summarize(day0 = unique(day0),
#             n = unique(n),
#             n0 = unique(n0)) %>% 
#   ungroup() %>% 
#   #filter(n0 > 0, n > 0) %>% 
#   mutate(log_r = log(n + 1) - log(n0 + 1))
# 
# g_nr <- df_n %>% 
#   mutate(key = case_when(light == 50 & treatment == "pc" & replicate == "c" ~ "Outlier",
#                          light == 100 & treatment == "pc" & replicate == "c" ~ "Outlier",
#                          light == 100 & treatment == "cp" & replicate %in% c("b", "c") ~ "Outlier",
#                          light == 100 & treatment == "cp" & replicate %in% c("a") ~ "Excluded",
#                          light == 200 & treatment == "cp" & replicate %in% c("a", "b", "c") ~ "Excluded",
#                          TRUE ~ as.character("Included"))) %>% 
#   ggplot(aes(x = n0,
#              y = log_r,
#              color = key)) +
#   geom_point() +
#   # facet_wrap(facets = ~treatment + light,
#   #            ncol = 4, nrow = 2,
#   #            scales = "free") +
#   facet_grid(rows = vars(light),
#              cols = vars(treatment, replicate),
#              scales = "free") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank())
