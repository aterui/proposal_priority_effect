
// State-space Gompertz model
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  // data:
  DATA_VECTOR(y);
  DATA_VECTOR(x1);
  DATA_VECTOR(x2);
  
  // parameters:
  PARAMETER(b0); // population growth rate parameter
  PARAMETER(log_b1); // frequency dependence
  PARAMETER(log_b2); // frequency dependence
  PARAMETER(log_sigma); // log(sigma)
  
  // procedures: (transformed parameters)
  Type b1 = exp(log_b1);
  Type b2 = exp(log_b2);
  Type sigma = exp(log_sigma);
  
  // reports on transformed parameters:
  ADREPORT(b1);
  ADREPORT(b2);
  ADREPORT(sigma);
  
  int n = y.size(); // get time series length
  
  Type nll = 0.0; // initialize negative log likelihood
  
  // model:
  for(int i = 1; i < n; i++){
    nll -= dnorm(y[i], b0 - b1 * x1[i] - b2 * x2[i], sigma, true);
  }
  
  return nll;
}
