
# set model ---------------------------------------------------------------

library(TMB)
compile("code/ricker_tmb.cpp")
dyn.load(dynlib("code/ricker_tmb"))


# simulated data ----------------------------------------------------------

sim_ricker <- function(N = 100,
                       seed = 123,
                       sigma_proc = 0.2,
                       sigma_obs = 0.5,
                       a = 1.5,
                       b = -0.01,
                       y1 = 4) {
  
  set.seed(seed)
  ytrue <- numeric(length = N)
  ytrue[1] <- y1
  proc_error <- rnorm(N, mean = 0, sd = sigma_proc)
  for(i in 2:N) {
    ytrue[i] <- ytrue[i - 1] + 
      a + 
      b * exp(ytrue[i - 1]) + 
      proc_error[i - 1]
  }
  x <- seq_len(N)
  yobs <- rnorm(N, mean = ytrue, sd = sigma_obs)
  y <- rpois(N, lambda = exp(yobs))
  
  return(y) 
}

y <- sim_ricker(N = 100, b = -0.005, seed = 1)
data <- list(y = y)

parameters <- list(r = 1,
                   b = 0.5,
                   log_sigma_proc = -1,
                   log_sigma_obs = -1,
                   u_obs = rep(mean(log(y)), length(y)),
                   u = rep(mean(log(y)), length(y)))
obj <- MakeADFun(data, parameters, random = c("u", "u_obs"), DLL = "ricker_tmb")
obj$hessian <- FALSE
opt <- nlminb(start = obj$par, obj = obj$fn, gr = obj$gr)
(rep <- sdreport(obj))

summary(rep, "fixed")
