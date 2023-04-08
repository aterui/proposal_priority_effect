pacman::p_load(tidyverse,
               foreach,
               cdyns)


x <- foreach(i = 1:100,
             .combine = bind_rows) %do% {
  
  print(i)
  if (i > 50) alpha <- 0 else alpha <- 1
  
  list_dyn <- cdynsim(n_timestep = 20,
                      n_species = 6,
                      alpha = alpha,
                      immigration = 1,
                      sd_immigration = 1,
                      model = "bh")
  
  df0 <- list_dyn$df_dyn %>% 
    mutate(y1 = rpois(nrow(.), lambda = density)) %>% 
    group_by(timestep) %>% 
    mutate(y_hat = sum(y1),
           f1 = y1 / y_hat) %>% 
    group_by(species) %>% 
    mutate(y0 = lag(y1),
           f0 = lag(f1)) %>% 
    ungroup() %>% 
    mutate(lambda = ifelse(y0 == 0,
                           NA,
                           y1 / y0),
           inv_f0 = ifelse(f0 == 0,
                           NA,
                           1 / f0),
           f_sp = factor(species)) %>% 
    drop_na(y0, lambda, inv_f0)
  
  # df0 %>%
  #   ggplot(aes(y = lambda,
  #              x = f0,
  #              color = f_sp)) +
  #   geom_point() +
  #   geom_smooth(method = "lm", se = F, linewidth = 0.1)# +
  #   #facet_wrap(facets =~ species)

  fit0 <- lm(lambda ~ inv_f0, df0)
  fit1 <- lm(lambda ~ f0 + inv_f0, df0)

  # fit0 <- glmmTMB::glmmTMB(lambda ~ (1 | species) + inv_f0, df0)
  # fit1 <- glmmTMB::glmmTMB(lambda ~ (1 | species) + f0 + inv_f0, df0)
  
  delta <- AIC(fit0) - AIC(fit1)
  
  return(tibble(alpha, delta))
}


x %>% 
  ggplot(aes(x = delta,
             color = factor(alpha),
             fill = factor(alpha))) +
  geom_density(alpha = 0.1)
