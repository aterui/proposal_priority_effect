
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))
#compile("code/lm.cpp")

cl <- makeCluster(detectCores() - 4)
registerDoSNOW(cl)


# sim data ----------------------------------------------------------------

df_para <- expand.grid(nsp = c(5, 10, 20),
                       nt = c(20, 40, 80),
                       alpha = seq(0, 1, by = 0.2),
                       sigma_alpha = c(0, 0.0001, 0.1, 0.25),
                       sigma_proc = 0.1,
                       min_k = c(500, 1000),
                       max_k = c(500, 1000),
                       min_r = c(0.5, 1.5, 2.5),
                       max_r = c(0.5, 1.5, 2.5),
                       neutral = c(0, 1)) %>% 
  mutate(m = min_k * 0.1/100) %>% # immigration: 0.1% carrying cap
  as_tibble() %>% 
  filter(max_k == min_k,
         max_r >= min_r,
         !(neutral == 1 & alpha != 1),
         !(neutral == 0 & alpha == 1),
         !(neutral == 1 & sigma_alpha != 0),
         !(neutral == 0 & sigma_alpha == 0)) %>% 
  mutate(group = row_number())

# test --------------------------------------------------------------------

pb <- txtProgressBar(max = nrow(df_para), style = 3)
fun_progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = fun_progress)

n_rep <- 100

tic()
df_sim <- foreach(x = iterators::iter(df_para, by = "row"),
                  .combine = bind_rows, 
                  .packages = c("tidyverse",
                                "foreach"),
                  .options.snow = opts) %dopar% {
                    
                    df0 <- foreach(j = 1:n_rep,
                                   .combine = bind_rows) %do% {
                                     
                                     tryCatch({
                                       ## set parameter values
                                       nsp <- x$nsp
                                       nt <- x$nt
                                       alpha <- x$alpha
                                       
                                       v_r <- runif(nsp, x$min_r, x$max_r)
                                       v_k <- runif(nsp, x$min_k, x$min_k)
                                       
                                       if (x$neutral == 1) {
                                         
                                         alpha <- 1
                                         A <- matrix(1, nrow = nsp, ncol = nsp)
                                         v_r <- mean(c(x$min_r, x$max_r))
                                         v_k <- mean(c(x$min_k, x$max_k))
                                         
                                       } else {
                                         A <- matrix(rbeta(nsp^2,
                                                           shape1 = alpha / x$sigma_alpha,
                                                           shape2 = (1 - alpha) / x$sigma_alpha),
                                                     nrow = nsp,
                                                     ncol = nsp)
                                         
                                         diag(A) <- 1
                                       }
                                       
                                       ## run simulation; output df_b
                                       list_dyn <- cdyns::cdynsim(n_species = nsp,
                                                                  n_timestep = nt,
                                                                  r_type = "constant",
                                                                  r = v_r,
                                                                  int_type = "manual",
                                                                  alpha = A,
                                                                  k = v_k,
                                                                  seed = 10,
                                                                  sd_env = x$sigma_proc,
                                                                  model = "ricker",
                                                                  immigration = x$m)
                                       
                                       ## add observation error
                                       df0 <- list_dyn$df_dyn %>%
                                         mutate(x = rpois(nrow(.), lambda = density))
                                       
                                       ## total community abundance
                                       df_n <- df0 %>% 
                                         group_by(timestep) %>% 
                                         summarize(n0 = sum(x)) %>% 
                                         ungroup() %>% 
                                         filter(timestep != max(timestep)) %>% 
                                         mutate(timestep = timestep + 1)
                                       
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
                                         mutate(lambda = x / x0,
                                                lambda = ifelse(is.infinite(lambda),
                                                                NA, # remove 0 counts from log_r estimate
                                                                lambda)) %>% 
                                         drop_na(lambda) %>% 
                                         left_join(df_n, by = "timestep")
                                       
                                       ## TMB fitting
                                       dyn.load(TMB::dynlib("code/lm"))
                                       parameters <- list(b0 = 1,
                                                          log_b1 = log(0.5),
                                                          log_sigma = -1)
                                       
                                       df_coef <- lapply(1:n_distinct(df_r$species),
                                                         function(i) {
                                                           df_i <- filter(df_r, species == i)
                                                           data <- list(y = df_i$lambda,
                                                                        x1 = df_i$p0)
                                                           
                                                           obj <- TMB::MakeADFun(data, parameters, DLL = "lm")
                                                           obj$hessian <- FALSE
                                                           opt <- nlminb(start = obj$par,
                                                                         obj = obj$fn,
                                                                         gr = obj$gr)
                                                           
                                                           return(list(species = i,
                                                                       b0 = opt$par[1],
                                                                       log_b1 = opt$par[2],
                                                                       message = opt$message)
                                                           )
                                                         }) %>% 
                                         bind_rows()
                                       
                                       
                                       ## mean frequency
                                       df_p <- df0 %>% 
                                         group_by(timestep) %>% 
                                         mutate(total = sum(x),
                                                p = x / total) %>% 
                                         ungroup() %>% 
                                         group_by(species) %>% 
                                         summarize(log_p = mean(log(p)),
                                                   log_p = ifelse(is.infinite(log_p),
                                                                  NA,
                                                                  log_p))
                                       
                                       ## combine
                                       df_b <- df_coef %>% 
                                         left_join(df_p,
                                                   by = "species")
                                       
                                       ## return
                                       return(tibble(df_b, x, replicate = j))
                                     }, error = function(x) cat("ERROR :",conditionMessage(x), "\n"))
                                     
                                   }
                    
                    return(df0)
                  }
toc()

stopCluster(cl)

df_z <- df_sim %>% 
  group_by(group, replicate) %>% 
  do(c = coef(lm(log_b1 ~ log_p, data = .))[1],
     z = coef(lm(log_b1 ~ log_p, data = .))[2]) %>% 
  ungroup() %>% 
  mutate(across(.cols = where(is.list), .fns = unlist),
         z_dev = abs(z - (-2))) %>% # deviation from -2
  left_join(df_para, by = "group")

# export ------------------------------------------------------------------

saveRDS(df_sim, file = "output/data_sim_tmb.rds")
saveRDS(df_z, file = "output/data_exponent_tmb.rds")
