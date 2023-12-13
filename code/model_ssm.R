model {
  
  tau0 <- 0.1
  
  # prior -------------------------------------------------------------------
  for (i in 1:Nsp) {
    for (j in 1:Ntr) {
      log_x0[i, j] ~ dt(0, tau0, 6)
      log_x[1, i, j] <- log_x0[i, j]
      g[i, j] ~ dnorm(0, tau0)T(-10, 10)
      alpha[i, j] ~ dnorm(0, tau0)
      beta[i, j] ~ dnorm(0, tau0)
    }
  }
  
  for (k in 1:2) {
    tau[k] ~ dscaled.gamma(2.5, 6)
    sigma[k] <- sqrt(1 / tau[k])
  }

  # likelihood --------------------------------------------------------------
  for (n in 1:Nsample) {
    Y[n] ~ dpois(x_obs[Timestep[n], Species[n], Treat[n]])
  }
    
  for (i in 1:Nsp) { # n species
    for (j in 1:Ntr) { # treatment
      
      # observation process
      # note: t decodes the first day of both species present as 1  
      for (t in 1:Nt[j]) { # time series
        log(x_obs[t, i, j]) <- log_x_obs[t, i, j]
        log_x_obs[t, i, j] ~ dnorm(log_x[t, i, j], tau[1])
        log(x[t, i, j]) <- log_x[t, i, j]
      } #t
      
      # state process
      # note: t decodes the first day of both species present as 1  
      for (t in 1:(Nt[j] - 1)) {
        log_x[t + 1, i, j] ~ dnorm(log_mu_x[t, i, j], tau[2])
        log_mu_x[t, i, j] <- 
          log_x[t, i, j] + 
          #log(int[t, i, j]) +
          g[i, j] + 
          alpha[i, j] * x[t, i, j] + 
          beta[i, j] * sum(x[t, 1:Nsp, j])
      }
      
    }
  }
  
  # re-parameterize ---------------------------------------------------------
  
  for (j in 1:Ntr) {
    mu_a[j] <- mean(alpha[,j])
    mu_b[j] <- mean(beta[,j])
  }
  
}

data {
  
  for (n in Index_int) {
    int[Timestep[n], Species[n], Treat[n]] <- Int[n]
  }
  
}