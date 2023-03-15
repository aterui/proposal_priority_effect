
// State-space Gompertz model
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
// data:
DATA_VECTOR(y);

// parameters:
PARAMETER(r); // population growth rate parameter
PARAMETER(b); // density dependence parameter
PARAMETER(log_sigma_proc); // log(process SD)
PARAMETER(log_sigma_obs); // log(observation SD)
PARAMETER_VECTOR(u_obs); // unobserved state vector
PARAMETER_VECTOR(u); // unobserved state vector

// procedures: (transformed parameters)
Type sigma_proc = exp(log_sigma_proc);
Type sigma_obs = exp(log_sigma_obs);

// reports on transformed parameters:
ADREPORT(sigma_proc)
ADREPORT(sigma_obs)
  
int n = y.size(); // get time series length

Type nll = 0.0; // initialize negative log likelihood

// process model:
for(int i = 1; i < n; i++){
  Type m = u[i - 1] + r + b * exp(u[i - 1]); // Ricker
  nll -= dnorm(u[i], m, sigma_proc, true);
}

// observation model:
for(int i = 0; i < n; i++){
  nll -= dnorm(u_obs[i], u[i], sigma_obs, true);
  nll -= dpois(y[i], exp(u_obs[i]), true);
}

return nll;
}