
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/function.R")
const <- 1

# data --------------------------------------------------------------------

## base data ####
df0 <- read_csv("data_raw/data_ojima2022.csv") %>% 
  mutate(across(.cols = starts_with("sp"),
                .fns = function(x) str_replace_all(x, pattern = "\\.", "_") %>% 
                  str_remove_all(pattern = "\\s")),
         pair = paste0(sp1, "_", sp2),
         x1 = total_od * frac_sp1,
         x2 = total_od * (1 - frac_sp1))

## log_r ####
# df_r0 - base dataframe for log_r
# df_r - keep if sample size for each species were > 3

df_r0 <- df0 %>% 
  group_by(treatment, pair) %>% 
  mutate(x1_0 = lag(x1),
         x2_0 = lag(x2),
         log_r_x1 = ifelse(treatment == "simultaneous",
                           log(x1 + const) - log(x1_0 + const),
                           NA),
         log_r_x2 = ifelse(treatment == "simultaneous",
                           log(x2 + const) - log(x2_0 + const),
                           NA)) %>% 
  ungroup() %>% 
  mutate(across(.cols = starts_with("log_r"),
                .fns = function(z) ifelse(is.nan(z)|is.infinite(z), NA, z))) %>% 
  drop_na(x1_0, x2_0) %>% 
  pivot_longer(cols = starts_with("log_r"),
               names_to = "type",
               values_to = "log_r") %>% 
  mutate(response = ifelse(str_detect(type, "_x1"), "x1", "x2")) %>% 
  select(-type) %>% 
  pivot_longer(cols = c("x1", "x2"),
               names_to = "predictor",
               values_to = "z0") %>% 
  filter(response == predictor) %>% 
  mutate(n0 = x1_0 + x2_0) %>% 
  drop_na(log_r)

# for sample size check
rm_pair <- df_r0 %>% 
  group_by(pair, response) %>% 
  tally() %>% 
  filter(n < 4) %>% 
  pull(pair) %>% 
  unique()

# remove small sample size data
df_r <- df_r0 %>% 
  filter(!(pair %in% rm_pair))

# community
df_n <- df0 %>% 
  filter(treatment == "simultaneous") %>% 
  group_by(pair) %>% 
  mutate(n = total_od,
         n0 = lag(total_od),
         log_r = log(n + const) - log(n0 + const)) %>% 
  drop_na(log_r) %>% 
  ungroup() %>% 
  filter(!(pair %in% rm_pair))

## strength of priority effects
df_gap <- df0 %>% 
  mutate(frac_sp2 = 1 - frac_sp1) %>% 
  filter(treatment != "simultaneous",
         step >= 4) %>% 
  group_by(treatment, pair) %>% 
  summarize(x1 = max(x1),
            x2 = max(x2)) %>% 
  pivot_wider(id_cols = pair,
              names_from = "treatment",
              values_from = c(x1, x2)) %>% 
  mutate(rr_x = abs(log(x1_sp1_first) - log(x1_sp2_first)),
         rr_y = abs(log(x2_sp2_first) - log(x2_sp1_first)),
         log_rr = 0.5 * (abs(rr_x) + abs(rr_y))) %>% 
  ungroup()

# analysis ----------------------------------------------------------------

v_pair <- unique(df_r$pair)

df_delta0 <- foreach(i = seq_len(length(v_pair)),
                     .combine = bind_rows) %do% {
                      
                      # observation
                      df_sub <- df_r %>%
                        filter(pair == v_pair[i])
                      
                      m <- MASS::rlm(log_r ~ n0 + z0 * response,
                                     data = df_sub,
                                     maxit = 2000)
                      
                      cout <- tibble(pair = v_pair[i],
                                     delta = coef(m)[2])
                      
                      # null
                      df_sub_n <- df_n %>%
                        filter(pair == v_pair[i])
                      
                      m0 <- MASS::rlm(log_r ~ n0,
                                      data = df_sub_n,
                                      maxit = 2000)
                      
                      r_hat <- coef(m0)[1]
                      a_hat <- coef(m0)[2]
                      se_a <- sqrt(diag(vcov(m0)))[2]
                      
                      p_neg <- fGarch::pstd(0, mean = a_hat, sd = se_a, nu = nrow(df_sub_n))
                      
                      if (a_hat < 0) {
                        set.seed(i)
                        delta0 <- null(n_species = n_distinct(df_sub$response),
                                       n_timestep = max(df_sub$step),
                                       r_hat = r_hat,
                                       a_hat = a_hat,
                                       seed = 50,
                                       sd_env = sd(resid(m0) * m0$w),
                                       nsim = 100,
                                       const = const)
                        
                        cout <- cout %>%
                          mutate(nsim = length(delta0),
                                 n = sum(cout$delta < delta0),
                                 p = mean(cout$delta < delta0),
                                 n_sample = nrow(df_sub_n),
                                 p_neg = p_neg)
                        
                      } else {
                        
                        cout <- cout %>%
                          mutate(nsim = NA,
                                 n = NA,
                                 p = NA,
                                 n_sample = nrow(df_sub_n),
                                 p_neg = p_neg)
                        
                      }
                      
                      return(cout)
                    }

df_delta <- df_delta0 %>% 
  left_join(df_gap,
            by = "pair")

saveRDS(df_delta,
        file = "output/simulation_ojima.rds")

# figure ------------------------------------------------------------------

df_delta %>% 
  ggplot(aes(x = p,
             y = log_rr)) +
  geom_point()
