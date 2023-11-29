
# setup -------------------------------------------------------------------

source("code/library.R")
source("code/function.R")

df_param <- expand.grid(n_timestep = c(10, 30),
                        n_species = c(5, 15),
                        r = 1,
                        sd_r = c(0, 0.25),
                        a0 = c(0.01, 0.02),
                        lambda_a1 = c(0, 0.5, 1, 1.5),
                        lambda_sd_a1 = c(0, 0.25),
                        n_rep = 100) %>% 
  as_tibble() %>% 
  mutate(a1 = a0 * lambda_a1,
         sd_a1 = a1 * lambda_sd_a1)


# simulation --------------------------------------------------------------

cl <- makeCluster(detectCores() - 4)
registerDoParallel(cl)

df_sim <- foreach(x = iterators::iter(df_param, by = "row"),
                  .packages = c("tidyverse", "foreach", "cdyns"),
                  .combine = bind_rows) %dopar% {
                    
                    df_out <- foreach(k = iterators::icount(25),
                                      .combine = bind_rows) %do% {
                                        
                                        cout <-sim(n_timestep = x$n_timestep,
                                                   n_species = x$n_species,
                                                   r = x$r,
                                                   sd_r = x$sd_r,
                                                   a0 = x$a0,
                                                   a1 = x$a1,
                                                   sd_a1 = x$sd_a1,
                                                   n_rep = x$n_rep)
                                        
                                        return(tibble(replicate = k, x, p = cout$p, b = cout$b))
                                      }
                    
                    return(df_out)
                  }

stopCluster(cl); gc()

saveRDS(df_sim, "output/simulation.rds")
