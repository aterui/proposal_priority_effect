
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/function.R")


# get lambda --------------------------------------------------------------

df_param <- expand.grid(n_species = c(2, 5, 10),
                        nsim = 100,
                        x0 = c(2, 8, 32),
                        a0 = 1 / c(2, 8, 32),
                        factor_h_x0 = seq(0, 1, by = 0.5),
                        factor_a1 = seq(0, 2, by = 0.2),
                        factor_h_a1 = seq(0, 0.5, by = 0.1)) %>%
  as_tibble() %>%
  mutate(a1 = round(a0 * factor_a1, 10),
         h_a1 = a1 * factor_h_a1,
         h_x0 = factor_h_x0)

df_out <- foreach(i = seq_len(nrow(df_param)),
                  .combine = bind_rows) %do% {
                    
                    print(i)
                    x <- df_param[i,]
                    
                    A <- with(x, {
                      if (h_a1 > 0) {
                        v_a1 <- runif(n = n_species^2,
                                      min = max(a1 - h_a1, 0),
                                      max = a1 + h_a1)
                      } else {
                        v_a1 <- rep(a1, n_species^2)
                      }
                      
                      ma <- matrix(v_a1, n_species, n_species)
                      diag(ma) <- a0
                      return(ma)
                    })
                    
                    x0_prime <- with(x,
                                     runif(n = n_species,
                                           min = x0 - h_x0,
                                           max = x0 + h_x0)
                    )
                    
                    x0_prime[x0_prime < 0] <- 0
                    
                    lambda <- with(x,
                                   stability(n_species = n_species,
                                             x0 = x0_prime,
                                             A = A))
                    
                    mu_r <- ifelse(is.na(lambda),
                                   NA,
                                   mean(attr(lambda, "R")))
                    
                    max_r <- ifelse(is.na(lambda),
                                   NA,
                                   max(attr(lambda, "R")))
                    
                    sd_r <- ifelse(is.na(lambda),
                                   NA,
                                   sd(attr(lambda, "R")))
                                   
                    xout <- x %>% 
                      mutate(lambda = lambda,
                             mu_r = mu_r,
                             max_r = max_r,
                             sd_r = sd_r,
                             n_persist = length(x0_prime[x0_prime > 0]))
                  }


set.seed(10)
df_set0 <- df_out %>% 
  filter(mu_r < 2) %>% 
  mutate(stability = ifelse(lambda < 1,
                            "insensitive",
                            "sensitive")) %>% 
  drop_na(stability) %>% 
  group_by(n_species, stability) %>% 
  sample_n(100) %>% 
  ungroup()

df_set <- bind_rows(mutate(df_set0, n_timestep = 10),
                    mutate(df_set0, n_timestep = 30)) %>% 
  relocate(n_timestep)

# figure ------------------------------------------------------------------

df_set %>% 
  ggplot(aes(x = lambda)) +
  geom_histogram() + 
  facet_grid(cols = vars(n_timestep),
             rows = vars(n_species))

# export ------------------------------------------------------------------

saveRDS(df_out, file = "data_fmt/param_set_src.rds")
saveRDS(df_set, file = "data_fmt/param_set.rds")
