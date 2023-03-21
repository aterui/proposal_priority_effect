
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))

cl <- makeCluster(detectCores() - 4)
registerDoSNOW(cl)


# sim data ----------------------------------------------------------------

df_para <- expand.grid(alpha = seq(0, 1, by = 0.25),
                       nsp = c(5, 10, 20),
                       nt = c(20, 50),
                       sigma_alpha = c(0, 0.1, 0.25),
                       sigma_proc = 0.05,
                       min_k = c(500, 1000),
                       max_k = c(500, 1000),
                       min_r = c(0.5, 1.5, 2.5),
                       max_r = c(0.5, 1.5, 2.5)) %>% 
  as_tibble() %>% 
  filter(max_k >= min_k,
         max_r >= min_r) %>% 
  mutate(group = row_number())


# test --------------------------------------------------------------------

pb <- txtProgressBar(max = nrow(df_para), style = 3)
fun_progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = fun_progress)

n_rep <- 100
n_maxit <- 1000

tic()
df_sim <- foreach(x = iterators::iter(df_para, by = "row"),
                  .combine = bind_rows, 
                  .packages = c("tidyverse",
                                "foreach"),
                  .options.snow = opts) %dopar% {
                    
                    df0 <- foreach(j = 1:n_rep,
                                   .combine = bind_rows) %do% {
                                     
                                     ## set parameter values
                                     nsp <- x$nsp
                                     nt <- x$nt
                                     
                                     v_r <- runif(nsp, x$min_r, x$max_r)
                                     k <- runif(nsp, x$min_k, x$min_k)
                                     
                                     A <- matrix(runif(nsp^2,
                                                       min = max(x$alpha - x$sigma_alpha, 0),
                                                       max = x$alpha + x$sigma_alpha),
                                                 nrow = nsp,
                                                 ncol = nsp)
                                     
                                     diag(A) <- 1
                                     
                                     if(x$alpha == 1) {
                                       A[,] <- 1
                                       k <- mean(k)
                                     }
                                     
                                     ## run simulation; output df_b
                                     list_dyn <- cdyns::cdynsim(n_species = nsp,
                                                                n_timestep = nt,
                                                                r_type = "constant",
                                                                r = v_r,
                                                                int_type = "manual",
                                                                alpha = A,
                                                                k = k,
                                                                seed = 100,
                                                                sd_env = x$sigma_proc,
                                                                model = "ricker",
                                                                immigration = 10)
                                     
                                     ## add observation error
                                     df0 <- list_dyn$df_dyn %>%
                                       mutate(x = rpois(nrow(.), lambda = density))
                                     
                                     ## calc frequency dependence
                                     df_r <- df0 %>% 
                                       group_by(timestep) %>% 
                                       mutate(total = sum(x),
                                              p = x / total) %>% 
                                       ungroup() %>% 
                                       group_by(species) %>% 
                                       mutate(x0 = lag(x),
                                              p0 = lag(p)) %>% 
                                       drop_na(x0, p0) %>% 
                                       ungroup() %>% 
                                       mutate(log_r = log(x) - log(x0),
                                              log_r = ifelse(is.infinite(log_r),
                                                             NA, # remove 0 counts from log_r estimate
                                                             log_r)) %>% 
                                       drop_na(log_r)
                                     
                                     ## TMB fitting
                                     dyn.load(TMB::dynlib("code/lm"))
                                     parameters <- list(b0 = 1,
                                                        log_b1 = log(0.5),
                                                        log_sigma = -1)
                                     
                                     df_coef <- lapply(1:n_distinct(df_r$species),
                                            function(i) {
                                              df_i <- filter(df_r, species == i)
                                              data <- list(y = df_i$log_r, x = df_i$p0)
                                              
                                              obj <- TMB::MakeADFun(data, parameters, DLL = "lm")
                                              obj$hessian <- FALSE
                                              opt <- nlminb(start = obj$par,
                                                            obj = obj$fn,
                                                            gr = obj$gr)
                                              
                                              return(list(species = i,
                                                          b0 = opt$par[1],
                                                          b1 = exp(opt$par[2])))
                                            }) %>% 
                                       bind_rows()
                                     
                                     
                                     ## mean frequency
                                     df_p <- df0 %>% 
                                       mutate(total = sum(x)) %>% 
                                       group_by(species) %>% 
                                       summarize(sum_x = sum(x),
                                                 total = unique(total)) %>% 
                                       ungroup() %>% 
                                       mutate(p = sum_x / total) %>% 
                                       dplyr::select(species, p)
                                     
                                     ## combine
                                     df_b <- df_coef %>% 
                                       left_join(df_p,
                                                 by = "species")
                                     
                                     ## return
                                     return(tibble(df_b, x, replicate = j))
                                     
                                   }
                    
                    return(df0)
                  }
toc()

stopCluster(cl)

df_z <- df_sim %>% 
  group_by(group, replicate) %>% 
  do(const = coef(MASS::rlm(log(b1) ~ log(p),
                            data = .,
                            maxit = n_maxit))[1],
     z = coef(MASS::rlm(log(b1) ~ log(p),
                        data = .,
                        maxit = n_maxit))[2]) %>% 
  ungroup() %>% 
  mutate(across(.cols = where(is.list), .fns = unlist),
         z_dev = abs(z - (-2))) %>% # deviation from -2
  left_join(df_para, by = "group")

# export ------------------------------------------------------------------

saveRDS(df_sim, file = "output/data_sim.rds")
saveRDS(df_z, file = "output/data_exponent.rds")
