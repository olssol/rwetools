//
//  POWER PRIOR ONLY FOR CONTINUOUS DATA
//  Power parameter A is fixed
//

data {
  //existing data
  int<lower = 1>  N0;
  real            YBAR0;
  real<lower = 0> SD0;

  //target borrowing
  real<lower = 0> A;
}

transformed data {
  real<lower = 0> sn0;
  sn0 = SD0/sqrt(N0*1.0);
}

parameters {
  real          theta;
}

model {
  //prior
  theta ~ normal(0, 1000);

  target +=  normal_lpdf(YBAR0 | theta, sn0) * A;
}
