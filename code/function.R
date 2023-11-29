
sim <- function(n_timestep,
                n_species,
                r,
                sd_r,
                a0,
                a1,
                sd_a1,
                sd_env = 0.1,
                n_rep = 100) {
  
  source("code/library.R")
  
  v_r <- runif(n_species, min = r - sd_r, max = r + sd_r)
  
  v_a1 <- abs(rnorm(n_species^2, mean = a1, sd = sd_a1))
  ma <- matrix(v_a1, n_species, n_species)
  diag(ma) <- a0
  
  list_dyn <- cdynsim(n_timestep = n_timestep,
                      n_species = n_species,
                      r = v_r,
                      alpha = ma,
                      int_type = "manual",
                      alpha_scale = "unscaled",
                      sd_env = sd_env)
  
  ## fitting
  df_dyn <- list_dyn$df_dyn %>% 
    group_by(timestep) %>% 
    mutate(nt1 = sum(density)) %>% 
    ungroup() %>% 
    group_by(species) %>% 
    mutate(n0 = lag(density),
           n1 = density,
           log_r = log(n1 / n0),
           nt0 = lag(nt1)) %>% 
    ungroup() %>% 
    mutate(species = factor(species))
  
  df_c <- df_dyn %>% 
    distinct(timestep, nt0, nt1) %>% 
    mutate(log_r = log(nt1 / nt0))
  
  fit <- lm(log_r ~ n0 + nt0, data = df_dyn)
  b <- coef(fit)[3]
  
  fit_c <- lm(log_r ~ nt0,
              data = df_c)
  
  sd_env <- sd(resid(fit_c))
  
  ## simulated null distribution
  r_hat <- coef(fit_c)[1]
  ma_hat <- matrix(-coef(fit_c)[2], n_species, n_species)
  
  v_beta <- foreach(i = iterators::icount(n_rep), .combine = c) %do% {
    
    list_sim <- cdynsim(n_timestep = n_timestep,
                        n_species = n_species,
                        r = r_hat,
                        alpha = ma_hat,
                        int_type = "manual",
                        alpha_scale = "unscaled",
                        sd_env = sd_env)
    
    df_sim <- list_sim$df_dyn %>% 
      group_by(timestep) %>% 
      mutate(nt1 = sum(density)) %>% 
      ungroup() %>% 
      group_by(species) %>% 
      mutate(n0 = lag(density),
             n1 = density,
             log_r = log(n1 / n0),
             nt0 = lag(nt1)) %>% 
      ungroup() %>% 
      mutate(species = factor(species))
    
    fit_sim <- lm(log_r ~ n0 + nt0, data = df_sim)
    beta0 <- coef(fit_sim)[3]
    
    return(beta0)
  }
  
  output <- list(p = mean(b < v_beta),
                 b = b,
                 b_null = v_beta)
  
  attr(output, "A") <- ma  
  attr(output, "R") <- v_r  
  
  return(output)
}


partial <- function(r, a, i, x0, model = "ricker") {
  
  ## vectorized parameters and variables
  v_a <- paste0("a[", seq_len(length(a)), "]")
  v_x <- paste0("x[", seq_len(length(x0)), "]")
  arg <- paste(c("r", "x", "a"), collapse = ", ")
  
  ## linear combination
  lcm <- paste(v_a, "*", v_x)
  
  if (model == "ricker") {
    
    fm <- c("r", lcm)
    m <- paste0("x[", i, "]", " * ",
                "exp(",
                paste(fm, collapse = " - "),
                ")")
    
    f_text <- parse(text = paste("f <- function(", arg, ") {",
                                 m,
                                 "}"))
    
    eval(f_text)
  }
  
  if (model == "bh") {
    
    m <- paste0("x[", i, "]", " * ", "exp(r)",
                " * ",
                "(1 + ",
                paste(lcm, collapse = " + "),
                ") ** -1")
    
    f_text <- parse(text = paste("f <- function(", arg, ") {",
                                 m,
                                 "}"))
    
    eval(f_text)
  }
    
  return(pracma::jacobian(f, x0 = x0, r = r, a = a))
}

stability <- function(n_species, R, A, model = "ricker") {
  
  # check input -------------------------------------------------------------
  
  if (any(unique(dim(A)) != n_species))
    stop("dimension mismatch in A")
  
  if (length(R) != n_species)
    stop("dimension mismatch in A")
  

  # get maximum absolute eigen value ----------------------------------------

  if (det(A) == 0) {
    
    return(NA)
    
  } else {
    
    if (model == "ricker") {
      x0 <- drop(solve(A) %*% R)
      J <- t(sapply(seq_len(n_species),
                    function(i) {
                      partial(r = R[i],
                              a = A[i, ],
                              x0 = x0,
                              i = i,
                              model = model)
                    }))
      
      lambda <- eigen(J)
      max_lambda <- max(abs(lambda$values))
    }
    
    if (model == "bh") {
      x0 <- drop(solve(A) %*% (exp(R) - 1))
      J <- t(sapply(seq_len(n_species),
                    function(i) {
                      partial(r = R[i],
                              a = A[i, ],
                              x0 = x0,
                              i = i,
                              model = model)
                    }))
      
      lambda <- eigen(J)
      max_lambda <- max(abs(lambda$values))
    }
    
    attr(max_lambda, "J") <- J
    
    return(max_lambda)
  }
  
}
