//
//  TYPICAL POWER PRIOR FOR CONTINUOUS DATA
//  Power parameter A is fixed
//

data {
  real<lower=0, upper=1> A;
  int  N0;
  int  N1;
  int<lower=0, upper=1> Y0[N0];
  int<lower=0, upper=1> Y1[N1];
}

parameters {
  real<lower=0, upper=1> theta;
}

model {
  //prior

  //likelihood
  if (N0 > 0) {
    target +=  bernoulli_lpmf(Y0 | theta)  * A;    
  }

  target +=  bernoulli_lpmf(Y1 | theta);
}
