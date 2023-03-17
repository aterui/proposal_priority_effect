
list_dyn <- cdyns::cdynsim(n_species = nsp,
                           n_timestep = nt,
                           r_type = "constant",
                           r = v_r,
                           int_type = "manual",
                           alpha = A,
                           k = k,
                           seed = 100,
                           sd_env = 0.1,
                           model = "ricker",
                           immigration = k/100)

df0 <- list_dyn$df_dyn %>%
  mutate(x = rpois(nrow(.), lambda = density))

df_b <- df0 %>% group_by(timestep) %>% 
  mutate(total = sum(x),
         p = x / total) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  mutate(x0 = lag(x),
         p0 = lag(p)) %>% 
  drop_na(x0, p0) %>% 
  mutate(log_r = log(x) - log(x0)) %>% 
  do(b0 = coef(MASS::rlm(log_r ~ p0, data = .))[1],
     b1 = coef(MASS::rlm(log_r ~ p0, data = .))[2]) %>% 
  ungroup() %>% 
  mutate(across(.cols = where(is.list), .fns = unlist),
         b1 = ifelse(b1 > 0, 0, b1))

df_p <- df0 %>% 
  mutate(total = sum(x)) %>% 
  group_by(species) %>% 
  summarize(sum_x = sum(x),
            total = unique(total)) %>% 
  ungroup() %>% 
  mutate(p = sum_x / total) %>% 
  dplyr::select(species, p)

df_b <- df_b %>% 
  left_join(df_p,
            by = "species")
