
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/function.R")


# get lambda --------------------------------------------------------------

df_param <- expand.grid(n_species = c(5, 15),
                        n_timestep = c(10, 30),
                        nsim = 100,
                        x0 = c(1, 5, 10),
                        factor_h_x0 = c(0, 0.5, 1),
                        a0 = c(0.01, 0.1),
                        factor_a1 = seq(0, 2, by = 0.25),
                        factor_h_a1 = seq(0, 0.5, by = 0.1)) %>% 
  as_tibble() %>%  
  mutate(a1 = round(a0 * factor_a1, 10),
         h_a1 = a1 * factor_h_a1,
         h_x0 = factor_h_x0)

df_out <- foreach(i = seq_len(nrow(df_param)),
                  .combine = bind_rows) %do% {
                    
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
                    
                    sd_r <- ifelse(is.na(lambda),
                                   NA,
                                   sd(attr(lambda, "R")))
                                   
                    xout <- x %>% 
                      mutate(lambda = lambda,
                             r = mu_r,
                             sd_r = sd_r,
                             n_persist = length(x0_prime[x0_prime > 0]))
                  }


# figure ------------------------------------------------------------------

df_out %>% 
  ggplot(aes(x = factor_a1,
             y = lambda,
             color = factor(a0))) +
  geom_point() +
  facet_grid(rows = vars(factor_h_x0),
             cols = vars(factor_h_a1, x0, n_species))


# export ------------------------------------------------------------------

df_param_set <- df_out %>% 
  filter(r < 2) %>% 
  mutate(stability = ifelse(lambda > 1, 1, 0)) %>% 
  group_by(n_species, n_timestep, stability) %>% 
  sample_n(150) %>% 
  ungroup()

saveRDS(df_param_set, file = "data_fmt/param_set.rds")
