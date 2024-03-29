
source("code/library.R")

# simulation --------------------------------------------------------------

sim <- function(n_timestep,
                n_species,
                x0,
                h_x0,
                a0,
                a1,
                h_a1,
                sd_env = 0.05,
                nsim = 100,
                const = 0) {
  
  if (h_a1 > 0) {
    v_a1 <- runif(n = n_species^2,
                  min = max(a1 - h_a1, 0),
                  max = a1 + h_a1)
  } else {
    v_a1 <- rep(a1, n_species^2)
  }
  
  ma <- matrix(v_a1, n_species, n_species)
  diag(ma) <- a0
  
  v_x0 <- runif(n_species,
                min = x0 - h_x0,
                max = x0 + h_x0)
  
  v_x0[v_x0 < 0] <- 0
  
  v_r <- drop(ma %*% v_x0)
  
  list_dyn <- suppressMessages(
    cdynsim(n_timestep = n_timestep,
            n_species = n_species,
            r = v_r,
            alpha = ma,
            int_type = "manual",
            alpha_scale = "unscaled",
            sd_env = sd_env,
            immigration = const)
  )
  
  ## data for linear fitting
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
    mutate(species = factor(species)) %>% 
    drop_na(log_r)
  
  df_c <- df_dyn %>% 
    distinct(timestep, nt0, nt1) %>% 
    mutate(log_r = log(nt1 / nt0)) %>% 
    drop_na(log_r)
  
  if (n_distinct(df_dyn$species) != n_species) {
    # NOTE - community grow infinitely; exclude
    output <- list(p = NA,
                   b = NA,
                   b_null = NA,
                   n = NA)
    
  } else {
    
    ## parameter observed
    fit <- lm(log_r ~ n0 + nt0 + species, data = df_dyn)
    b <- coef(fit)[3]
    
    ## parameter for a null model
    fit_c <- lm(log_r ~ nt0, data = df_c)
    sd_env <- sd(resid(fit_c))
    r_hat <- coef(fit_c)[1]
    ma_hat <- matrix(-coef(fit_c)[2], n_species, n_species)
    
    v_beta <- foreach(i = iterators::icount(nsim),
                      .combine = c) %do% {
      
      list_sim <- suppressMessages(
        cdynsim(n_timestep = n_timestep,
                n_species = n_species,
                r = r_hat,
                alpha = ma_hat,
                int_type = "manual",
                alpha_scale = "unscaled",
                sd_env = sd_env,
                immigration = const)
      )
      
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
        mutate(species = factor(species)) %>%
        drop_na(log_r)
      
      if (n_distinct(df_sim$species) != n_species) {
        beta0 <- NA
      } else {
        fit_sim <- lm(log_r ~ n0 + nt0 + species, data = df_sim)
        beta0 <- coef(fit_sim)[3]
      }
      
      return(beta0)
    }
    
    v_beta <- na.omit(v_beta)
    names(v_beta) <- NULL
    p <- mean(b < v_beta)
    output <- list(p = ifelse(is.nan(p), NA, p),
                   b = b,
                   b_null = v_beta,
                   n = length(v_beta))
    
  }
  
  attr(output, "A") <- ma  
  attr(output, "R") <- v_r  
  
  return(output)
}


# data simulation ---------------------------------------------------------

simdata <- function(n_timestep,
                    n_species,
                    x0,
                    h_x0,
                    a0,
                    a1,
                    h_a1,
                    sd_env = 0.05,
                    const = 0) {
  
  if (h_a1 > 0) {
    v_a1 <- runif(n = n_species^2,
                  min = max(a1 - h_a1, 0),
                  max = a1 + h_a1)
  } else {
    v_a1 <- rep(a1, n_species^2)
  }
  
  ma <- matrix(v_a1, n_species, n_species)
  diag(ma) <- a0
  
  v_x0 <- runif(n_species,
                min = x0 - h_x0,
                max = x0 + h_x0)
  
  v_x0[v_x0 < 0] <- 0
  
  v_r <- drop(ma %*% v_x0)
  
  list_dyn <- suppressMessages(
    cdynsim(n_timestep = n_timestep,
            n_species = n_species,
            r = v_r,
            alpha = ma,
            int_type = "manual",
            alpha_scale = "unscaled",
            sd_env = sd_env,
            immigration = const)
  )
  
  ## data for linear fitting
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
    mutate(species = factor(species)) %>% 
    drop_na(log_r)
  
  df_c <- df_dyn %>% 
    distinct(timestep, nt0, nt1) %>% 
    mutate(log_r = log(nt1 / nt0)) %>% 
    drop_na(log_r)
  
  return(list(data_i = df_dyn,
              data_c = df_c,
              A = ma,
              r = v_r))  
}

# partial -----------------------------------------------------------------

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


# stability ---------------------------------------------------------------

stability <- function(n_species, R, x0, A, model = "ricker") {
  
  # check input -------------------------------------------------------------
  
  if (any(unique(dim(A)) != n_species))
    stop("dimension mismatch in A")
  
  if (!missing(R)) {
    if (length(R) != n_species)
      stop("dimension mismatch in A or R")
  }
  
  if (!missing(x0)) {
    if (length(x0) != n_species)
      stop("dimension mismatch in A or x0")
  }
  
  # get maximum absolute eigen value ----------------------------------------

  if (det(A) == 0) {
    
    return(NA)
    
  } else {
    
    if (model == "ricker") {
      
      if (missing(R) & !missing(x0)) {
        R <- drop(A %*% x0)
      }
      
      if (!missing(R) & missing(x0)) {
        x0 <- drop(solve(A) %*% R)
      }
      
      # check negative equilibrium
      if (any(x0 < 0)) {
        # species going extinct
        e <- which(x0 < 0)

        # remove species that goes extinct
        A_prime <- A[-e, -e]
        R_prime <- R[-e]
        # re-calculate equilibrium
        x0_prime <- solve(A_prime) %*% R_prime

        # update equilibrium value
        x0[-e] <- x0_prime
        x0[e] <- 0
      }
      
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
      
      if (missing(R) & !missing(x0)) {
        R <- drop(log(1 + A %*% x0))
      }
      
      if (!missing(R) & missing(x0)) {
        x0 <- drop(solve(A) %*% (exp(R) - 1))
      }

      # check negative equilibrium
      if (any(x0 < 0)) {
        # species going extinct
        e <- which(x0 < 0)
        
        # remove species that goes exinct
        A_prime <- A[-e, -e]
        R_prime <- R[-e]
        # re-calculate equilibrium
        x0_prime <- solve(A_prime) %*% R_prime
        
        # update equilibrium value
        x0[-e] <- x0_prime
        x0[e] <- 0
      }
      
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
    attr(max_lambda, "R") <- R
    attr(max_lambda, "x0") <- x0
    
    return(max_lambda)
  }
  
}


# null model function -----------------------------------------------------

null <- function(n_species, n_timestep, r_hat, a_hat, sd_env, 
                 seed = 50,
                 nsim = 100,
                 const = 1) {
  
  ma_hat <- matrix(-a_hat, n_species, n_species)
  
  v_beta <- foreach(i = iterators::icount(nsim),
                    .combine = c) %do% {
    
    list_sim <- suppressMessages(
      cdyns::cdynsim(n_timestep = n_timestep,
                     n_species = n_species,
                     #n_warmup = 0,
                     #n_burnin = 0,
                     r = r_hat,
                     alpha = ma_hat,
                     int_type = "manual",
                     alpha_scale = "unscaled",
                     sd_env = sd_env,
                     seed = seed,
                     immigration = 0)
    )
    
    df_sim <- list_sim$df_dyn %>% 
      group_by(timestep) %>% 
      mutate(n = sum(density)) %>% 
      ungroup() %>% 
      group_by(species) %>% 
      mutate(z0 = lag(density),
             z = density,
             log_r = log(z + const) - log(z0 + const),
             n0 = lag(n)) %>% 
      ungroup() %>% 
      mutate(species = factor(species)) %>% 
      drop_na(log_r)
    
    if (n_distinct(df_sim$species) == n_species) {
      fit_sim <- lm(log_r ~ n0 + z0 * species, data = df_sim)
      beta0 <- coef(fit_sim)[2]
    } else {
      beta0 <- NA
    }
    
    return(beta0)
  }
  
  return(v_beta)
}
