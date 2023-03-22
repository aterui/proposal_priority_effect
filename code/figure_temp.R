
# setup -------------------------------------------------------------------
rm(list = ls())
source(here::here("code/library.R"))

df_b <- readRDS(here::here("output/data_exponent.rds"))

# plot --------------------------------------------------------------------

df_b %>%
  filter(sigma_alpha == 0.25,
         nt == 20,
         nsp == 10,
         min_k == 250,
         max_k == min_k) %>% 
  ggplot(aes(x = factor(alpha),
             y = log(abs(z + 2)))) +
  #geom_violin(draw_quantiles = 0.5) +
  geom_jitter(alpha = 0.1) +
  facet_grid(rows = vars(max_r),
             cols = vars(min_r),
             scales = "free") +
  #geom_hline(yintercept = 0) +
  theme_bw() +
  theme(panel.grid = element_blank())
  
