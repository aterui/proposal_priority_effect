#+ include = F
# setup -------------------------------------------------------------------
rm(list = ls())
source(here::here("code/library.R"))

df_b <- readRDS(here::here("output/data_exponent.rds"))


# plot --------------------------------------------------------------------

#' 1. No variation in alpha, r, k
#+ echo = F, message = F, warning = F

df_b %>%
  filter(sigma_alpha == 0,
         min_r == max_r) %>% 
  ggplot(aes(x = factor(alpha),
             y = z_dev)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  facet_wrap(facets = ~ min_r + nsp,
             scales = "free",
             labeller = label_both)

#' 2. No variation in alpha; Variation in r, k
#+ echo = F, message = F, warning = F

df_b %>%
  filter(sigma_alpha == 0,
         min_r == 0.5,
         max_r == 2.5) %>% 
  ggplot(aes(x = factor(alpha),
             y = z_dev)) +
  geom_boxplot() +
  #geom_jitter(alpha = 0.2) +
  facet_wrap(facets = ~ min_r + nsp,
             scales = "free",
             labeller = label_both)

#' 3. Variation in alpha; No variation in r, k
#+ echo = F, message = F, warning = F, fig.height = 8

df_b %>%
  filter(sigma_alpha != 0,
         min_r == max_r) %>% 
  ggplot(aes(x = factor(alpha),
             y = z_dev)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  facet_wrap(facets = ~ min_r + nsp,
             scales = "free",
             labeller = label_both)


#' 4. Variation in alpha, r, k; SD alpha = 0.1
#+ echo = F, message = F, warning = F

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
