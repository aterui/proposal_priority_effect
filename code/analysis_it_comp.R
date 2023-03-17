
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))


# sim data ----------------------------------------------------------------

df_para <- expand.grid(alpha = seq(0, 1, by = 0.25),
                       nsp = c(5, 10, 20),
                       nt = 20,
                       sigma_alpha = c(0, 0.1, 0.25),
                       min_k = c(500, 1000),
                       max_k = c(1000, 1500),
                       min_r = c(0.5, 1.5, 2.5),
                       max_r = c(0.5, 1.5, 2.5)) %>% 
  as_tibble() %>% 
  filter(max_k >= min_k,
         max_r >= min_r) %>% 
  mutate(group = row_number())

# test --------------------------------------------------------------------

n_rep <- 1
df_para <- slice(df_para, 1:2)

df_b <- foreach(x = iterators::iter(df_para, by = "row"),
                .combine = bind_rows) %do% {
                  
                  df0 <- foreach(j = 1:n_rep,
                                 .combine = bind_rows) %do% {
                                   
                                   print(paste0("group:", x$group, ", replicate:", j))
                                   
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
                                     v_r <- rep(mean(v_r), nsp)
                                     k <- mean(k)
                                   }
                                   
                                   source(here::here("code/sim_data.R"))
                                   
                                   return(tibble(df_b, x))
                                   
                                 }
                  
                  return(df0)
                }

saveRDS(df_b, file = "output/sim_re.rds")
