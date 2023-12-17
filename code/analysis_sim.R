
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/function.R")

df_param <- readRDS("data_fmt/param_set.rds") %>% 
  as_tibble()
sim_run <- purrr::possibly(sim, otherwise = NULL)

# simulation --------------------------------------------------------------

cl <- makeCluster(parallel::detectCores() - 2)
registerDoSNOW(cl)

pb <- txtProgressBar(max = nrow(df_param), style = 3)
fun_progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = fun_progress)

tictoc::tic()
df_sim <- foreach(i = seq_len(nrow(df_param)),
                  .packages = c("tidyverse", "foreach", "cdyns"),
                  .combine = bind_rows,
                  .options.snow = opts) %dopar% {
                    
                    x <- df_param[i, ]
                    
                    df_out <- foreach(k = iterators::icount(5),
                                      .combine = bind_rows) %do% {
                                        
                                        seed <- 100 * i + k
                                        set.seed(seed)
                                        
                                        cout <-with(x,
                                                    sim_run(n_timestep = n_timestep,
                                                            n_species = n_species,
                                                            x0 = x0,
                                                            h_x0 = h_x0,
                                                            a0 = a0,
                                                            a1 = a1, 
                                                            h_a1 = h_a1,
                                                            nsim = nsim))
                                        
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
                                                      eigen_max = eigen_max,
                                                      seed = seed))
                                      }
                    
                    return(df_out)
                  }
tictoc::toc()

stopCluster(cl); gc()

saveRDS(df_sim, "output/simulation.rds")
