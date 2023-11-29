
# setup -------------------------------------------------------------------

source("code/library.R")
source("code/function.R")

df_param <- expand.grid(n_timestep = c(10, 30),
                        n_species = c(5, 15),
                        r = 1,
                        sd_r = c(0, 0.25),
                        a0 = c(0.01, 0.02),
                        factor_a1 = c(0, 0.5, 1, 1.5),
                        factor_sd_a1 = c(0, 0.25),
                        n_rep = 100) %>% 
  as_tibble() %>% 
  mutate(a1 = a0 * factor_a1,
         sd_a1 = a1 * factor_sd_a1)


# simulation --------------------------------------------------------------

cl <- makeCluster(detectCores() - 4)
registerDoParallel(cl)

tictoc::tic()
df_sim <- foreach(x = iterators::iter(df_param %>% slice(1), by = "row"),
                  .packages = c("tidyverse", "foreach", "cdyns"),
                  .combine = bind_rows) %dopar% {
                    
                    df_out <- foreach(k = iterators::icount(5),
                                      .combine = bind_rows) %do% {
                                        
                                        cout <-with(x,
                                                    sim(n_timestep = n_timestep,
                                                        n_species = n_species,
                                                        r = r,
                                                        sd_r = sd_r,
                                                        a0 = a0,
                                                        a1 = a1, 
                                                        sd_a1 = sd_a1,
                                                        n_rep = n_rep))
                                        
                                        A <- attr(cout, "A")
                                        R <- attr(cout, "R")
                                        
                                        eigen_max <- stability(n_species = x$n_species,
                                                               R = R,
                                                               A = A)
                                        
                                        return(tibble(replicate = k,
                                                      x,
                                                      p = cout$p,
                                                      b = cout$b,
                                                      eigen_max = eigen_max))
                                      }
                    
                    return(df_out)
                  }
tictoc::toc()

stopCluster(cl); gc()

saveRDS(df_sim, "output/simulation.rds")
