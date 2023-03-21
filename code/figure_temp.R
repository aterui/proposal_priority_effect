
# setup -------------------------------------------------------------------
rm(list = ls())
source(here::here("code/library.R"))

df_b <- readRDS(here::here("output/data_exponent.rds"))


# plot --------------------------------------------------------------------

## 1. No variation in alpha, r, k
df_b %>%
  filter(sigma_alpha == 0,
         max_r == min_r,
         nt == 20,
         max_k > min_k) %>% 
  ggplot(aes(x = factor(alpha),
             y = z)) +
  geom_violin(draw_quantiles = 0.5) +
  facet_grid(rows = vars(nsp),
             cols = vars(min_r),
             scales = "free") +
  geom_hline(yintercept = -2) +
  theme_bw() +
  theme(panel.grid = element_blank())
  

df_b %>%
  filter(sigma_alpha == 0.1,
         nt == 20) %>% 
  ggplot(aes(x = sqrt(z_dev),
             y = alpha)) +
  geom_point(alpha = 0.05) +
  # geom_smooth(method = "glm",
  #             method.args = list(family = "binomial")) +
  facet_grid(rows = vars(nsp),
             cols = vars(min_r),
             scales = "free") +
  theme_bw()


## 2. No variation in alpha; Variation in r, k
df_b %>%
  filter(sigma_alpha == 0,
         min_r == 0.5,
         max_r == 2.5) %>% 
  ggplot(aes(x = factor(alpha),
             y = z_dev)) +
  geom_violin(draw_quantiles = 0.5) +
  facet_wrap(facets = ~ min_r + nsp,
             scales = "free",
             labeller = label_both)

## 3. Variation in alpha; No variation in r, k
df_b %>%
  filter(sigma_alpha == 0.1,
         min_r < max_r) %>% 
  ggplot(aes(x = factor(alpha),
             y = z_dev)) +
  geom_violin(draw_quantiles = 0.5) +
  facet_wrap(facets = ~ nsp,
             scales = "free",
             labeller = label_both)


## 4. Variation in alpha, r, k; SD alpha = 0.1
df_b %>%
  filter(min_r == 0.5, max_r == 1.5,
         min_k == 500, max_k == 1000) %>%
  ggplot(aes(x = factor(a),
             y = abs(value + 2))) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  facet_grid(rows = vars(sigma,
                         min_r, max_r,
                         min_k, max_k),
             cols = vars(nsp),
             scales = "free",
             labeller = label_both)
