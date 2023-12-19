
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/function.R")
const <- 1

# data --------------------------------------------------------------------

df0 <- read_csv("data_raw/data_hsu2021.csv",
                guess_max = 10000) %>% 
  rename_with(.fn = str_to_lower) %>% 
  dplyr::select(-growth_rate, -priority_effects, -chl_a)  %>% 
  rename(day = experimental_day,
         light = light_level) %>% 
  mutate(across(.cols = where(is.character),
                .fns = function(x) {
                  str_to_lower(x) %>% 
                    str_replace_all(pattern = "\\s",
                                    replacement = "_")
                }),
         light = case_when(light == "0" ~ 0,
                           light == "50" ~ 50,
                           light == "100" ~ 100,
                           light == "200" ~ 200)) %>% 
  filter(treatment %in% c("cp", "pc")) %>% 
  pivot_wider(id_cols = c(day, replicate, treatment, light),
              names_from = species,
              values_from = cells) %>% 
  rename(x = p_bur,
         y = colp) %>% 
  mutate(n = x + y) %>% 
  drop_na(n) %>% 
  group_by(replicate, treatment, light) %>% 
  mutate(key = ifelse(n > 0, 1, 0),
         last_day = ifelse(any(key == 0),
                           min(day[which(key == 0)]),
                           max(day))
  ) %>% 
  filter(day < last_day) %>% 
  ungroup()

df_r0 <- df0 %>% 
  group_by(replicate, treatment, light) %>% 
  arrange(day, by_group = TRUE) %>% 
  mutate(x0 = lag(x),
         y0 = lag(y),
         n0 = lag(n),
         log_r_x = log(x + const) - log(x0 + const),
         log_r_y = log(y + const) - log(y0 + const),
         log_r_n = log(x + const) - log(x0 + const),
         day0 = lag(day),
         interval = day - day0) %>%
  mutate(across(.cols = starts_with("log_r"),
                .fns = function(z) ifelse(is.nan(z) | is.infinite(z), NA, z))) %>% 
  ungroup() %>% 
  relocate(day, day0, interval) %>% 
  drop_na(day0) %>% 
  pivot_longer(cols = starts_with("log_r"),
               names_to = "response",
               values_to = "log_r") %>% 
  pivot_longer(cols = c("x0", "y0"),
               names_to = "predictor",
               values_to = "z0") %>% 
  mutate(predictor = str_remove_all(predictor, "0")) %>% 
  separate(col = response,
           into = c("null1",
                    "null2",
                    "species"),
           sep = "_") %>% 
  dplyr::select(-starts_with("null")) %>% 
  filter(species == predictor)

df_n <- df_r0 %>%
  filter(species == "x") %>% 
  mutate(log_r = log(n + const) - log(n0 + const)) %>% 
  select(-c(x, y, species, predictor, z0))

df_r <- df_r0 %>% 
  drop_na(log_r)

saveRDS(df_r, file = "data_fmt/data_hsu_fmt.rds")

# analysis ----------------------------------------------------------------

df_i <- df_r %>% 
  distinct(replicate, treatment, light)

df_delta <- foreach(i = seq_len(nrow(df_i)),
                    .combine = bind_rows) %do% {
                      
                      # observation
                      df_sub <- df_r %>%
                        filter(replicate == df_i$replicate[i],
                               treatment == df_i$treatment[i],
                               light == df_i$light[i])
                      
                      m <- MASS::rlm(log_r ~ n0 + z0 * species + offset(log(interval)),
                                     data = df_sub,
                                     maxit = 2000)
                      
                      cout <- with(df_i,
                                   tibble(replicate = replicate[i],
                                          treatment = treatment[i],
                                          light = light[i],
                                          delta = coef(m)[2]))

                      # null
                      df_sub_n <- df_n %>%
                        filter(replicate == df_i$replicate[i],
                               treatment == df_i$treatment[i],
                               light == df_i$light[i])

                      m0 <- MASS::rlm(log_r ~ n0 + offset(log(interval)),
                                      data = df_sub_n,
                                      maxit = 2000)
                      
                      r_hat <- coef(m0)[1]
                      a_hat <- coef(m0)[2]
                      se_a <- sqrt(diag(vcov(m0)))[2]
                      
                      p_neg <- fGarch::pstd(0, mean = a_hat, sd = se_a, nu = nrow(df_sub_n))
                      
                      if (a_hat < 0) {
                        set.seed(i)
                        delta0 <- null(n_species = n_distinct(df_sub$species),
                                       n_timestep = n_distinct(df_sub$day),
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

# strength of priority effects --------------------------------------------

# df_x: calculate the strength of priority effects as log_rr
#       x and y are densities of protists
df_x <- df0 %>% 
  group_by(replicate, treatment, light) %>% 
  summarize(x = max(x),
            y = max(y)) %>% 
  pivot_wider(id_cols = c("replicate", "light"),
              names_from = "treatment",
              values_from = c("x", "y")) %>% 
  ungroup() %>% 
  mutate(rr_x = abs(log(x_pc) - log(x_cp)),
         rr_y = abs(log(y_pc) - log(y_cp)),
         log_rr = 0.5 * (rr_x + rr_y)) %>% 
  arrange(light)

# merge with proposed statistic
df_delta <- df_delta %>% 
  left_join(df_x, by = c("replicate", "light"))

saveRDS(df_delta,
        file = "output/simulation_hsu.rds")

# figure ------------------------------------------------------------------

df_r %>%
  filter(species == "y") %>% 
  ggplot(aes(x = n0,
             y = log_r,
             color = factor(species))) +
  geom_point() +
  facet_grid(rows = vars(treatment, replicate),
             cols = vars(light),
             scales = "free")
