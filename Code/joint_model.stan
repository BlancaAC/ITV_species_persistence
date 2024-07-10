functions {

}
data {
  // lengths of different pieces of data
  int<lower=0> N_POI;
  int<lower=0> N_BIN;
  int<lower=1> N_POL;
  
  // response variables
  array[N_POI] int SEEDS_POI;
  array[N_BIN] int SUCCESS_BIN;
  array[N_BIN] int TRIAL_BIN;
  
  // predictor variables and covariates
  matrix[N_POI, N_POL] X_POI;
  matrix[N_BIN, N_POL] X_BIN;

  // to directly sample the prior
  int prior_only;
}
transformed data {

}
parameters {
  // uncontrained variable used to create postive constrained betas
  vector[N_POL] taus;
}
transformed parameters {
  // postive constrained betas
  vector<lower=0>[N_POL] betas;
  betas = log1p_exp(taus);

  // observation-wise predicted number of seeds per flower
  vector<lower=0>[N_POI] seeds_per_flower;
  seeds_per_flower = X_POI * betas;

  // observation-wise predicted probability of producing fruit
  vector<lower=0,upper=1>[N_BIN] probability;
  probability = -expm1(-X_BIN * betas);
}
model {
  // prior for unconstrained parameters tau
  taus ~ normal(0,2);

  if(!prior_only){
    // likelihood for POISSON observations
    target += poisson_lpmf(SEEDS_POI | seeds_per_flower);

    // likelihood for BERNOULLI observations
    target += binomial_lpmf(SUCCESS_BIN | TRIAL_BIN, probability);
  }
}
generated quantities {
  // calculate and save the log-likelihood
  real ll = 0;
  ll += poisson_lpmf(SEEDS_POI | seeds_per_flower);
  ll += binomial_lpmf(SUCCESS_BIN | TRIAL_BIN, probability);
}
