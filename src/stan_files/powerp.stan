//
//  TYPICAL POWER PRIOR FOR CONTINUOUS DATA
//  Power parameter A is fixed
//

data {
  real<lower=0, upper=1> A;
  int  N0;
  int  N1;
  real Y0[N0];
  real Y1[N1];
}

parameters {
  real          theta;
  real<lower=0> tau0;
  real<lower=0> tau1;
}

model {
  //prior
  theta ~ normal(0, 1000);
  tau0  ~ cauchy(0, 2.5);
  tau1  ~ cauchy(0, 2.5);

  //likelihood
  if (N0 > 0) {
    target +=  normal_lpdf(Y0 | theta, tau0) * A;
  }

  target +=  normal_lpdf(Y1 | theta, tau1);
}
