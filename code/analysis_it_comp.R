
# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))


# sim data ----------------------------------------------------------------

df_para <- expand.grid(a = seq(0, 1, by = 0.25),
                       nsp = c(5, 10, 20),
                       nt = 20,
                       sigma = c(0, 0.1, 0.25),
                       min_k = c(500, 1000),
                       max_k = c(1000, 1500),
                       min_r = c(0.5, 1.5, 2.5),
                       max_r = c(0.5, 1.5, 2.5)) %>% 
  as_tibble() %>% 
  filter(max_k >= min_k,
         max_r >= min_r)

# test --------------------------------------------------------------------

n_rep <- 50

df_b <- foreach(i = 1:nrow(df_para),
                .combine = bind_rows) %do% {
                  
                  x <- df_para[i, ]
                  nsp <- x$nsp
                  nt <- x$nt
                  
                  df0 <- foreach(j = 1:n_rep,
                                 .combine = bind_rows) %do% {
                                   
                                   print(paste0((i - 1) * n_rep + j, "/", n_rep * nrow(df_para)))
                                                                      
                                   v_r <- runif(nsp, x$min_r, x$max_r)
                                   k <- runif(nsp, x$min_k, x$min_k)
                                   
                                   A <- matrix(abs(rnorm(nsp^2,
                                                         mean = x$a,
                                                         sd = x$sigma)),
                                               nsp,
                                               nsp)
                                   
                                   diag(A) <- 1
                                   
                                   if(x$a == 1) {
                                     A[,] <- 1
                                     v_r <- rep(mean(v_r), nsp)
                                     k <- mean(k)
                                   }
                                   
                                   source(here::here("code/sim_data.R"))
                                   
                                   df_null <- df2 %>%
                                     filter(sp1 == sp2) %>% 
                                     group_by(t) %>%
                                     mutate(xt = sum(x_j),
                                            xt0 = sum(x0_j)) %>% 
                                     ungroup() %>% 
                                     mutate(p_x = x_i / xt,
                                            p_x0 = x0_i / xt0,
                                            log_pr = log(p_x) - log(p_x0))
                                   
                                   df_p <- list_dyn$df_species %>% 
                                     dplyr::select(species, mean_density) %>% 
                                     mutate(p = mean_density / sum(mean_density)) %>% 
                                     rename(sp1 = species)
                                   
                                   df_lm <- df_null %>%
                                     group_by(sp1) %>%
                                     do(lm = lm(log_pr ~ p_x0, .) %>% coef()) %>%
                                     mutate(b0 = lm[1],
                                            b1 = lm[2]) %>%
                                     dplyr::select(-lm) %>%
                                     left_join(df_p,
                                               by = "sp1")
                                   
                                   fit <- lm(log(abs(b1)) ~ log(p), df_lm) %>% 
                                     summary()
                                   
                                   return(
                                          tibble(x,
                                                 r = cor(df_lm$b1,
                                                         df_lm$p,
                                                         method = "spearman"),
                                                 rsq = fit$r.squared,
                                                 z = coef(fit)[2, 1])
                                          )
                                   
                                 }
                  
                  return(df0)
                }

saveRDS(df_b, file = "output/sim_re.rds")

df_b %>%
  pivot_longer(cols = c(r, rsq, z)) %>%
  ggplot(aes(x = factor(a),
             y = value)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_jitter(alpha = 0.2) +
  facet_wrap(facets = ~ name,
             scales = "free")
 
# df_null %>% 
#   ggplot(aes(x = p_x0,
#              y = log_pr,
#              color = factor(sp1))) +
#   geom_point() +
#   geom_smooth(method = "lm")
