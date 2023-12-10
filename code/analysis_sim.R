
# setup -------------------------------------------------------------------

source("code/library.R")
source("code/function.R")

df_param <- expand.grid(n_timestep = c(10, 30),
                        n_species = c(5, 10, 20),
                        r = 1,
                        sd_r = c(0, 0.1),
                        a0 = c(0.01, 0.05),
                        sd_a1 = c(0, 0.1),
                        factor_a1 = seq(0, 1.5, by = 0.25),
                        n_rep = 100) %>% 
  as_tibble() %>% 
  mutate(a1 = round(a0 * factor_a1, 10))

sim_run <- purrr::possibly(sim, otherwise = NULL)

# simulation --------------------------------------------------------------

cl <- makeCluster(parallel::detectCores() - 2)
registerDoSNOW(cl)

pb <- txtProgressBar(max = nrow(df_param), style = 3)
fun_progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = fun_progress)

tictoc::tic()
df_sim <- foreach(x = iterators::iter(df_param, by = "row"),
                  .packages = c("tidyverse", "foreach", "cdyns"),
                  .combine = bind_rows,
                  .options.snow = opts) %dopar% {
                    
                    df_out <- foreach(k = iterators::icount(3),
                                      .combine = bind_rows) %do% {
                                        
                                        cout <-with(x,
                                                    sim_run(n_timestep = n_timestep,
                                                            n_species = n_species,
                                                            r = r,
                                                            sd_r = sd_r,
                                                            a0 = a0,
                                                            a1 = a1, 
                                                            sd_a1 = sd_a1,
                                                            n_rep = n_rep))
                                        
                                        if (!is.null(cout)) {
                                          eigen_max <- stability(n_species = x$n_species,
                                                                 R = attr(cout, "R"),
                                                                 A = attr(cout, "A"))
                                        } else {
                                          cout <- NULL
                                          cout$p <- NA
                                          cout$b <- NA
                                          cout$n <- NA
                                        }
                                        
                                        return(tibble(replicate = k,
                                                      x,
                                                      p = cout$p,
                                                      b = cout$b,
                                                      n_sim = cout$n,
                                                      eigen_max = eigen_max))
                                      }
                    
                    return(df_out)
                  }
tictoc::toc()

stopCluster(cl); gc()

saveRDS(df_sim, "output/simulation.rds")
