
# setup -------------------------------------------------------------------
rm(list = ls())
source(here::here("code/library.R"))

df_alpha <- readRDS(here::here("output/data_sim.rds")) %>% 
  group_by(group, replicate) %>% 
  summarize(alpha = unique(alpha))

df_z <- readRDS(here::here("output/data_exponent.rds")) %>% 
  left_join(df_alpha, by = c("group", "replicate"))

# plot --------------------------------------------------------------------

df_z %>%
  filter(sigma_alpha == 0,
         nt == 20,
         nsp == 10,
         min_k == 100,
         max_k == min_k) %>% 
  ggplot(aes(x = alpha,
             y = abs(z + 2))) +
  geom_jitter(alpha = 0.1) +
  facet_grid(rows = vars(max_r),
             cols = vars(min_r),
             scales = "free") +
  geom_smooth(method = "loess",
              color = grey(0.1)) +
  geom_hline(yintercept = 1) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(panel.grid = element_blank())
  
