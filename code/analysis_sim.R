
# setup -------------------------------------------------------------------

source("code/library.R")
source("code/function.R")

df_param <- expand.grid(n_timestep = c(10, 30),
                        n_species = c(5, 15),
                        x0 = 1,
                        h_x0 = seq(0, 1, by = 0.5),
                        a0 = c(0.01, 0.05),
                        factor_a1 = seq(0, 1.5, by = 0.25),
                        factor_h_a1 = c(0.01, 0.25),
                        nsim = 100) %>% 
  as_tibble() %>% 
  mutate(a1 = round(a0 * factor_a1, 10),
         h_a1 = a1 * factor_h_a1) # rounded to avoid float point issue

sim_run <- purrr::possibly(sim, otherwise = NULL)

# simulation --------------------------------------------------------------

cl <- makeCluster(parallel::detectCores() - 2)
registerDoSNOW(cl)

pb <- txtProgressBar(max = nrow(df_param), style = 3)
fun_progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = fun_progress)

tictoc::tic()
df_sim <- foreach(i = 1:2,#seq_len(nrow(df_param)),
                  .packages = c("tidyverse", "foreach", "cdyns"),
                  .combine = bind_rows,
                  .options.snow = opts) %dopar% {
                    
                    x <- df_param[i, ]
                    
                    df_out <- foreach(k = iterators::icount(25),
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
