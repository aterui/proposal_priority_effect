
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))


# parameter ---------------------------------------------------------------

df_para <- expand.grid(alpha = seq(0, 1, by = 0.25),
                       nsp = c(5, 10, 20),
                       nt = 20,
                       sigma_alpha = c(0, 0.1, 0.25),
                       sigma_obs = 0.25,
                       sigma_proc = 0.1,
                       min_k = c(500, 1000),
                       max_k = c(1000, 1500),
                       min_r = c(0.5, 1.5, 2.5),
                       max_r = c(0.5, 1.5, 2.5)) %>% 
  as_tibble() %>% 
  filter(max_k >= min_k,
         max_r >= min_r) %>% 
  mutate(param_id = row_number())


# generate simulated data -------------------------------------------------

n_rep <- 2
x <- df_para[1,]

tic()
# df_simdata <- foreach(x = iterators::iter(df_para, by = "row"),
#                       .combine = bind_rows) %do% {
#                         
df0 <- foreach(j = 1:n_rep,
               .combine = bind_rows) %do% {
                 
                 ## parameter setup
                 v_r <- runif(x$nsp, x$min_r, x$max_r)
                 k <- runif(x$nsp, x$min_k, x$min_k)
                 
                 A <- matrix(abs(rnorm(x$nsp^2,
                                       mean = x$alpha,
                                       sd = x$sigma_alpha)),
                             x$nsp,
                             x$nsp)
                 
                 diag(A) <- 1
                 
                 if(x$alpha == 1) {
                   A[,] <- 1
                   v_r <- rep(mean(v_r), x$nsp)
                   k <- mean(k)
                 }
                 
                 ## simulate data
                 list_dyn <- cdyns::cdynsim(n_timestep = x$nt, 
                                            n_species = x$nsp,
                                            k = k,
                                            r = v_r,
                                            int_type = "manual",
                                            alpha = A,
                                            sd_env = x$sigma_proc)
                 
                 df_fit <- list_dyn$df_dyn %>% 
                   mutate(log_lambda = log(density) + rnorm(n = nrow(.),
                                                            mean = 0,
                                                            sd = x$sigma_obs),
                          count = rpois(n = nrow(.),
                                        lambda = exp(log_lambda))) %>% 
                   left_join(list_dyn$df_species,
                             by = c("species")) %>% 
                   arrange(species) %>% 
                   mutate(x)
                 
                 ## TMB: fit data to model
                 list_est <- lapply(seq_len(n_distinct(df_fit$species)),
                                    function(i) {
                                      
                                      y <- df_fit %>% 
                                        filter(species == i) %>% 
                                        pull(count)
                                      
                                      data <- list(y = y)
                                      
                                      params <- list(r = mean(v_r),
                                                     b = - mean(v_r) / mean(k),
                                                     log_sigma_proc = log(x$sigma_proc),
                                                     u = rep(mean(log(y)), length(y)))
                                      
                                      obj <- MakeADFun(data,
                                                       params,
                                                       random = c("u"),
                                                       DLL = "ricker_tmb")
                                      
                                      opt <- nlminb(start = obj$par, obj = obj$fn, gr = obj$gr)
                                      
                                      rep <- sdreport(obj)
                                      
                                    })
                 
                 n_t <- lapply(list_est,
                               function(x) summary(x, "random")[, 1]) %>% 
                   sapply(function(x) exp(x)) %>% 
                   apply(1, sum) %>% 
                   mean()
                 
                 b <- lapply(list_est,
                        function(x) summary(x, "fixed")[2, "Estimate"]) %>% 
                   unlist()
                 
                 df_p <- df_fit %>% 
                   group_by(species) %>% 
                   summarize(sum_i = sum(count)) %>% 
                   ungroup() %>% 
                   mutate(p = sum_i / sum(sum_i),
                          b = b,
                          z = n_t * b)
                 
                 
               }
#   
#   return(df0)
# }
toc()