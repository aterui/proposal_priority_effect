
// State-space Gompertz model
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  // data:
  DATA_VECTOR(y);
  DATA_VECTOR(x);
  
  // parameters:
  PARAMETER(b0); // population growth rate parameter
  PARAMETER(log_b1); // frequency dependence
  PARAMETER(log_sigma); // log(sigma)
  
  // procedures: (transformed parameters)
  Type b1 = exp(log_b1);
  Type sigma = exp(log_sigma);
  
  // reports on transformed parameters:
  ADREPORT(b1);
  ADREPORT(sigma);
  
  int n = y.size(); // get time series length
  
  Type nll = 0.0; // initialize negative log likelihood
  
  // model:
  for(int i = 1; i < n; i++){
    nll -= dnorm(y[i], b0 - b1 * x[i], sigma, true);
  }
  
  return nll;
}
