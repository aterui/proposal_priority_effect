
# setup -------------------------------------------------------------------
rm(list = ls())
source(here::here("code/library.R"))

#df_z <- readRDS(here::here("output/data_exponent_ssm.rds"))
df_z <- readRDS(here::here("output/data_exponent.rds")) %>% 
  mutate(sigma_alpha = factor(sigma_alpha))

df_frame <- df_z %>% 
  expand(distinct(., nsp, nt, min_r, max_r, min_k, max_k),
         sigma_alpha = factor(c(0.0001, 0.1, 0.25)))

df_null <- df_z %>% 
  filter(neutral == 1) %>% 
  drop_na(z_dev) %>% 
  group_by(nsp, nt, min_r, max_r, min_k, max_k) %>% 
  summarize(high = quantile(z_dev, 0.75),
            med = quantile(z_dev, 0.5),
            low = quantile(z_dev, 0.25)) %>% 
  ungroup() %>% 
  left_join(df_frame,
            by = c("nsp", "nt", "min_r", "max_r", "min_k", "max_k"),
            multiple = "all")
  

# plot --------------------------------------------------------------------

options(scipen=100)

## label
lab_min <- c(`0.5` = "r[min]==0.5",
             `1.5` = "r[min]==1.5",
             `2.5` = "r[min]==2.5")

lab_max <- c(`0.5` = "r[max]==0.5",
             `1.5` = "r[max]==1.5",
             `2.5` = "r[max]==2.5")

lab_sigma <- c(`1e-04` = "sigma[alpha]==10^{-4}",
               `0.1` = "sigma[alpha]==0.10",
               `0.25` = "sigma[alpha]==0.25")


## plot export
df_set <- distinct(df_z, nt, nsp, min_k)

foreach(x = iterators::iter(df_set, by = "row")) %do% {
                  
          df_z_plot <- df_z %>%
            filter(nt == x$nt,
                   nsp == x$nsp,
                   min_k == x$min_k)
          
          df_null_plot <- df_null %>%
            filter(nt == x$nt,
                   nsp == x$nsp,
                   min_k == x$min_k)
          
          ## plot
          g <- df_z_plot %>% 
            filter(neutral == 0) %>% 
            ggplot(aes(x = factor(alpha),
                       y = z_dev,
                       color = sigma_alpha,
                       fill = sigma_alpha)) +
            geom_hline(data = df_null_plot,
                       aes(yintercept = med),
                       linetype = "dashed",
                       color = grey(0.5)) +
            geom_hline(data = df_null_plot,
                       aes(yintercept = high),
                       linetype = "dotted",
                       color = grey(0.5)) +
            geom_hline(data = df_null_plot,
                       aes(yintercept = low),
                       linetype = "dotted",
                       color = grey(0.5)) +
            # geom_jitter(alpha = 0.25,
            #             size = 0.5) +
            geom_boxplot(alpha = 0.1,
                         linewidth = 0.1,
                         outlier.color = NA) +
            scale_y_continuous(trans = "log10") +
            facet_grid(rows = vars(sigma_alpha),
                       cols = vars(min_r, max_r),
                       #scales = "free",
                       labeller = labeller(min_r = as_labeller(lab_min,
                                                               label_parsed),
                                           max_r = as_labeller(lab_max,
                                                               label_parsed),
                                           sigma_alpha = as_labeller(lab_sigma,
                                                                     label_parsed))) +
            labs(x = expression(mu[alpha]),
                 y = expression(z[dev]),
                 color = expression(sigma[alpha]),
                 fill = expression(sigma[alpha])) +
            theme_bw() +
            theme(panel.grid = element_blank(),
                  strip.background = element_blank())
          
          ggsave(g,
                 filename = paste0("output/figure_lm_t_", x$nt, "_s_", x$nsp, "_k_", x$min_k, ".pdf"),
                 width = 12,
                 height = 8)
        }
