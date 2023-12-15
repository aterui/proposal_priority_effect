
# setup -------------------------------------------------------------------

source("code/function.R")


# get lambda --------------------------------------------------------------

df_param <- expand.grid(n_species = c(5, 15),
                        x0 = 1,
                        h_x0 = seq(0, 1, by = 0.5),
                        a0 = c(0.01, 0.05),
                        factor_a1 = seq(0, 1.5, by = 0.25),
                        factor_sd_a1 = c(0.01, 0.25),
                        nsim = 100) %>% 
  as_tibble() %>%  
  mutate(a1 = round(a0 * factor_a1, 10),
         h_a1 = a1 * factor_sd_a1)

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
                    
                    
                    lambda <- with(x,
                                   stability(n_species = n_species,
                                             x0 = runif(n_species,
                                                        max(0, x0 - h_x0),
                                                        x0 + h_x0),
                                             A = A))
                    
                    xout <- x %>% 
                      mutate(lambda = lambda,
                             r = mean(attr(lambda, "R")),
                             sd_r = sd(attr(lambda, "R")))
                  }


# figure ------------------------------------------------------------------

df_out %>% 
  ggplot(aes(x = factor_a1,
             y = lambda,
             color = factor(a0))) +
  geom_point() +
  facet_grid(rows = vars(h_x0),
             cols = vars(factor_sd_a1, x0, n_species))
