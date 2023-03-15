#+ include = F
# setup -------------------------------------------------------------------
rm(list = ls())
source(here::here("code/library.R"))

df_b <- readRDS(here::here("output/sim_re.rds"))


# plot --------------------------------------------------------------------

#' 1. No variation in alpha, r, k
#+ echo = F, message = F, warning = F

df_b %>%
  pivot_longer(cols = z) %>%
  filter(sigma == 0,
         min_r == max_r,
         min_k == max_k) %>% 
  ggplot(aes(x = factor(a),
             y = abs(value + 2))) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  facet_grid(rows = vars(min_r, min_k),
             cols = vars(nsp),
             scales = "free",
             labeller = label_both)

#' 2. No variation in alpha; Variation in r, k
#+ echo = F, message = F, warning = F

df_b %>%
  pivot_longer(cols = z) %>%
  filter(sigma == 0,
         min_r == 0.5, max_r == 1.5,
         min_k < max_k) %>%
  ggplot(aes(x = factor(a),
             y = abs(value + 2))) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  facet_grid(rows = vars(min_r, min_k, max_r, max_k),
             cols = vars(nsp),
             scales = "free",
             labeller = label_both)

#' 3. Variation in alpha; No variation in r, k
#+ echo = F, message = F, warning = F, fig.height = 8

df_b %>%
  pivot_longer(cols = z) %>%
  filter(sigma != 0,
         min_r == max_r,
         min_k == max_k) %>%
  ggplot(aes(x = factor(a),
             y = abs(value + 2))) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  facet_grid(rows = vars(min_r, min_k, sigma),
             cols = vars(nsp),
             scales = "free",
             labeller = label_both)


#' 4. Variation in alpha, r, k; SD alpha = 0.1
#+ echo = F, message = F, warning = F

df_b %>%
  pivot_longer(cols = z) %>%
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
