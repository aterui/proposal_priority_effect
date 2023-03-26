
// State-space model
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  // data:
  DATA_VECTOR(y); // abundance
  DATA_SCALAR(phi); // total community abundance
  int n = y.size(); // get time series length
  
  // parameters:
  PARAMETER(b0); // population growth rate parameter
  PARAMETER(b1); // log(frequency dependence parameter)
  PARAMETER(log_sigma_proc); // log(process SD)
  PARAMETER_VECTOR(u); // unobserved state vector; relative abundance
  Type nll = 0.0; // initialize negative log likelihood
  
  // procedures: (transformed parameters)
  Type sigma_proc = exp(log_sigma_proc);
  
  // reports on transformed parameters:
  ADREPORT(sigma_proc);
  
  // process model:
  for(int i = 1; i < n; i++){
    Type mu = u[i - 1] + b0 + b1 * u[i - 1]; // Ricker analogue
    nll -= dnorm(u[i], mu, sigma_proc, true);
  }
  
  // observation model:
  for(int i = 0; i < n; i++){
    Type lambda = exp(u[i]) * phi;
    nll -= dpois(y[i], lambda, true);
  }
  
  return nll;
}
