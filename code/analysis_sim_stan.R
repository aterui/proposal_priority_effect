
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/function.R")

df_param <- readRDS("data_fmt/param_set.rds") %>% 
  as_tibble()

cl <- makeCluster(parallel::detectCores() - 4)
registerDoSNOW(cl)

pb <- txtProgressBar(max = nrow(df_param), style = 3)
fun_progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = fun_progress)


# run simulation ----------------------------------------------------------

tictoc::tic()
df_sim <- foreach(i = seq_len(nrow(df_param)),
                  .packages = c("tidyverse", "foreach", "cdyns"),
                  .combine = bind_rows,
                  .options.snow = opts) %dopar% {
                    
                    x <- df_param[i, ]
                    
                    df_out <- foreach(k = iterators::icount(3),
                                      .combine = bind_rows) %do% {
                                        
                                        seed <- 100 * i + k
                                        set.seed(seed)
                                        
                                        # obtain simulated data
                                        cdata <-with(x,
                                                     simdata(n_timestep = n_timestep,
                                                             n_species = n_species,
                                                             x0 = x0,
                                                             h_x0 = h_x0,
                                                             a0 = a0,
                                                             a1 = a1, 
                                                             h_a1 = h_a1))
                                        
                                        # calculate eigen value
                                        eigen_max <- stability(n_species = x$n_species,
                                                               R = cdata$r,
                                                               A = cdata$A)
                                        
                                        # fit stan model to the data
                                        fit <- with(cdata,
                                                    rstanarm::stan_glmer(log_r ~ n0 + nt0 + (1 | species),
                                                                         data = data_i,
                                                                         chains = 3))
                                        
                                        # get scaled beta = beta / r by species
                                        mcmc <- as.matrix(fit)
                                        mu <- mcmc[, str_detect(colnames(mcmc), "^\\(Intercept\\)$")]
                                        eps <- mcmc[, str_detect(colnames(mcmc), "b\\[\\(Intercept\\).{1,}\\]")]
                                        b <- mcmc[, "nt0"]
                                        b_scale <- b / (mu + eps)
                                        
                                        # get scaled beta for community
                                        fit0 <- with(cdata,
                                                     rstanarm::stan_glm(log_r ~ nt0,
                                                                        data = data_c,
                                                                        chains = 3))
                                        
                                        b0_scale <- fit0$coefficients[2] / fit0$coefficients[1]
                                        
                                        # binary indicator s(b) vs. s(b0)
                                        # s(.) denotes scaled beta
                                        # repeat across mcmc samples
                                        n_exceed <- sapply(1:nrow(b_scale),
                                                           function(i) sum(b_scale[i, ] < b0_scale))
                                        
                                        p <- mean(n_exceed > 1)
                                        mu_n_exceed <- mean(n_exceed)
                                        
                                        return(tibble(replicate = k,
                                                      x,
                                                      n_exceed = mu_n_exceed,
                                                      p = p,
                                                      b = fit$coefficients["nt0"],
                                                      eigen_max = eigen_max,
                                                      seed = seed))
                                      }
                    
                    return(df_out)
                  }

tictoc::toc()


# export ------------------------------------------------------------------

stopCluster(cl); gc()
saveRDS(df_sim, "output/simulation_stan.rds")

