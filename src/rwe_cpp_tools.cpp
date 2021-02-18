#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::init]]
void my_package_init(DllInfo *dll) {
  // initialization code here
  R_useDynamicSymbols(dll, TRUE);
}

// expit to get probability
double g_ex(NumericVector beta, NumericVector x, double tol = 1e-6) {
  double ex;
  int k, n_beta = beta.size();

  ex = 0;
  for (k = 0; k < n_beta; k++) {
    ex += beta[k] * x(k);
  }
  ex = 1 / (1 + std::exp(-ex));

  if (ex < tol) {
    ex = tol;
  } else if (ex > 1 - tol) {
    ex = 1 - tol;
  }

  return(ex);
}


//' Get PS
//'
//' @export
// [[Rcpp::export]]
NumericVector c_get_ps(NumericVector beta, NumericMatrix mat_x, double tol = 1e-6) {
  int n_pat  = mat_x.nrow();
  NumericVector ps(n_pat);

  int i;

  for (i = 0; i < n_pat; i++) {
    ps[i] = g_ex(beta, mat_x(i, _), tol);
  }

  return(ps);
}


//' Test Rcpp function
//'
//'
//' @param test test parameter
//'
//' @export
// [[Rcpp::export]]
double crtTest(double test) {
  return pow(test,2);
}

//' Match by nearest neighbor
//'
//' Match subjects in group candidate with subject in the target group
//'
//' @param target  vector of data from the target group
//' @param candidates vector of data from the candidate group
//' @param ratio  1:ratio match
//'
// [[Rcpp::export]]
NumericVector cMatch(NumericVector target, NumericVector candidate, int ratio) {

  NumericVector dist(candidate.size());
  NumericVector inx(candidate.size());

  int i, j, start;

  // initialize index of candidates
  for (i = 0; i < inx.size(); i++) {
    inx[i] = i;
  }

  // find neighbors
  start = 0;
  for (i = 0; i < target.size(); i++) {

    // calculate distance
    for (j = start; j < inx.size(); j++) {
      dist[inx[j]] = fabs(target[i] - candidate[inx[j]]);
    }

    // sort distance
    std::sort(inx.begin() + start, inx.end(),
              [&](const int& a, const int& b) {
                return (dist[a] < dist[b]);
              });

    // keep ratio neighbors
    start += ratio;
  }

  // return
  return(inx);
}


//' Get Moment Constraints for Logistic Regression
void c_ps_gmm_g_logi(NumericMatrix rst,
                     IntegerVector grp,
                     NumericVector ps,
                     NumericMatrix mat_grp_x,
                     bool att = false) {

  // first column is group
  int n_beta = mat_grp_x.ncol() - 1;
  int n_pat  = mat_grp_x.nrow();

  double z, ex, w;
  int i, k;

  for (i = 0; i < n_pat; i++) {
    z  = grp[i];
    ex = ps[i];
    if (att) {
      w =  z * ex * (1 - ex) - (1 - z) * ex * ex;
    } else {
      w =  z * (1 - ex) - (1 - z) * ex;
    }

    for (k = 0; k < n_beta; k++) {
      rst(i, k) = w * mat_grp_x(i, k + 1);
    }
  }
}

//' Get Moment Constraints for balance in covariates
void c_ps_gmm_g_xbal(NumericMatrix rst,
                     IntegerVector grp,
                     NumericVector ps,
                     NumericMatrix mat_grp_x,
                     int offset = 0,
                     int moment = 1) {

  // first column is group
  int n_x   = mat_grp_x.ncol() - 2;
  int n_pat = mat_grp_x.nrow();

  NumericVector curx;

  double sw1 = 0, sw0 = 0, z, ex, w;
  int i, k;

  //sum of weights
  for (i = 0; i < n_pat; i++) {
    if (1 == grp[i]) {
      sw1 += 1 / ps[i];
    } else if (0 == grp[i]) {
      sw0 += 1 / (1 - ps[i]);
    }
  }

  //moments constraints
  for (i = 0; i < n_pat; i++) {
    curx  = mat_grp_x(i, _);
    z     = grp[i];
    ex    = ps[i];

    if (1 == z) {
      w = 1 / ex / sw1;
    } else if (0 == z) {
      w = - 1 / (1 - ex) / sw0;
    } else {
      w = 0;
    }

    for (k = 0; k < n_x; k++) {
      rst(i, offset + k) = w * pow(curx(k + 2), moment);
    }
  }
}


//' Get ALL Moment Constraints
//'
//' return N row k column matrix
//' Row:   subject
//' Column: Moment
//'
//'
// [[Rcpp::export]]
NumericMatrix c_ps_gmm_g(NumericVector beta, NumericMatrix mat_grp_x, bool att = false) {

  int n_pat  = mat_grp_x.nrow();
  int n_x    = mat_grp_x.ncol() - 2;
  int i;

  NumericVector ps(n_pat);
  IntegerVector mgrp(n_pat), grp(n_pat);
  NumericVector curx;

  NumericMatrix rst(n_pat, 5 * n_x + 1);
  //get integer group
  for (i = 0; i < n_pat; i++) {
    mgrp[i] = (int)(mat_grp_x(i, 0));
  }

  // get propensity score
  for (i = 0; i < n_pat; i++) {
    curx  = mat_grp_x(i, _);
    ps[i] = g_ex(beta, curx[Range(1, n_x + 2)]);
  }

  // logistic regression
  for (i = 0; i < n_pat; i++) {
    if (mgrp[i] == 1) {
      grp[i] = 1;
    } else {
      grp[i] = 0;
    }
  }
  c_ps_gmm_g_logi(rst, grp, ps, mat_grp_x, att);

  // x-bal between grp 1 vs. 2
  for (i = 0; i < n_pat; i++) {
    if (mgrp[i] == 1) {
      grp[i] = 1;
    } else if (mgrp[i] == 2){
      grp[i] = 0;
    } else {
      grp[i] = -1;
    }
  }
  c_ps_gmm_g_xbal(rst, grp, ps, mat_grp_x, n_x + 1);
  c_ps_gmm_g_xbal(rst, grp, ps, mat_grp_x, 2 * n_x + 1, 2);

  // x-bal between grp 1 vs. 3
  for (i = 0; i < n_pat; i++) {
    if (mgrp[i] == 1) {
      grp[i] = 1;
    } else if (mgrp[i] == 3) {
      grp[i] = 0;
    } else {
      grp[i] = -1;
    }
  }
  c_ps_gmm_g_xbal(rst, grp, ps, mat_grp_x, 3 * n_x + 1);
  c_ps_gmm_g_xbal(rst, grp, ps, mat_grp_x, 4 * n_x + 1, 2);

  // return
  return(rst);
}

//' Get Gbar and Sigma
//'
//' return a List
//'
//'
// [[Rcpp::export]]
NumericVector c_ps_gmm_gbar(NumericVector beta,
                            NumericMatrix mat_x,
                            NumericVector ps,
                            IntegerVector mgrp,
                            NumericVector mz,
                            NumericVector mr,
                            NumericVector n3) {

  int n_pat  = mat_x.nrow();
  int n_x    = mat_x.ncol() - 1;
  int dim_g  = 3 * n_x + 1;
  int i, j;
  double tmp, tmp2;

  NumericVector gbar(dim_g);

  //initialize
  std::fill(gbar.begin(), gbar.end(), 0);

  // g_bar: add all
  for (i = 0; i < n_pat; i++) {
    // logistic
    tmp = mz[i] * (1 - ps[i]) - (1 - mz[i]) * ps[i];
    for (j = 0; j < n_x + 1; j++) {
      gbar[j] += tmp * mat_x(i, j);
    }

    // grp = 1 vs. 2
    tmp  = mz[i] / ps[i] - (1 - mz[i]) * mr[i]       / (1 - ps[i]);
    tmp2 = mz[i] / ps[i] - (1 - mz[i]) * (1 - mr[i]) / (1 - ps[i]);

    for (j = 0; j < n_x; j++) {
      gbar[j + n_x     + 1] += tmp  * mat_x(i, j+1);
      gbar[j + 2 * n_x + 1] += tmp2 * mat_x(i, j+1);
    }
  }

  // g_bar: take average
  for (j = 0; j < n_x + 1; j++) {
    gbar[j] /= n_pat;
  }

  for (j = 0; j < n_x; j++) {
    gbar[j + n_x     + 1] /= (n3[0] + n3[1]);
    gbar[j + 2 * n_x + 1] /= (n3[0] + n3[2]);
  }

  // return
  return(gbar);
}

//' Get Gbar and Sigma
//'
//' return a List
//'
//' @export
// [[Rcpp::export]]
NumericMatrix c_ps_gmm_sigma(NumericVector beta,
                             NumericMatrix mat_x,
                             NumericVector ps,
                             IntegerVector mgrp,
                             NumericVector mz,
                             NumericVector mr,
                             NumericVector n3) {

  int n_pat  = mat_x.nrow();
  int n_x    = mat_x.ncol() - 1;
  int dim_g  = 3 * n_x + 1;
  int i, j, k;
  double tmp, tmp2;

  NumericMatrix sigma(dim_g, dim_g);
  NumericMatrix sigma11(n_x + 1, n_x + 1);
  NumericMatrix sigma22(n_x,     n_x);
  NumericMatrix sigma33(n_x,     n_x);
  NumericMatrix sigma12(n_x + 1, n_x);
  NumericMatrix sigma13(n_x + 1, n_x);
  NumericMatrix sigma23(n_x,     n_x);

  //initialize
  std::fill(sigma.begin(),   sigma.end(),   0);
  std::fill(sigma11.begin(), sigma11.end(), 0);
  std::fill(sigma22.begin(), sigma22.end(), 0);
  std::fill(sigma33.begin(), sigma33.end(), 0);
  std::fill(sigma12.begin(), sigma12.end(), 0);
  std::fill(sigma13.begin(), sigma13.end(), 0);
  std::fill(sigma23.begin(), sigma23.end(), 0);

  // Sigma: add all
  for (i = 0; i < n_pat; i++) {
    //sigma11 logistic
    tmp = ps[i] * (1 - ps[i]);
    for (j = 0; j < n_x + 1; j++) {
      for (k = j; k < n_x + 1; k++) {
        sigma11(j, k) += tmp * mat_x(i, j) * mat_x(i, k);
      }
    }

    //sigma12 sigma13
    for (j = 0; j < n_x + 1; j++) {
      for (k = 0; k < n_x; k++) {
        if (3 != mgrp[i]) {
          sigma12(j, k) += mat_x(i, j) * mat_x(i, k+1);
        }

        if (2 != mgrp[i]) {
          sigma13(j, k) += mat_x(i, j) * mat_x(i, k+1);
        }
      }
    }

    // sigma22 sigma33 sigma23
    tmp  = 1 / ps[i] / (1 - ps[i]);
    tmp2 = 1 / ps[i];
    for (j = 0; j < n_x; j++) {
      for (k = j; k < n_x; k++) {
        if (3 != mgrp[i]) {
          sigma22(j, k) += tmp * mat_x(i, j+1) * mat_x(i, k+1);
        }

        if (2 != mgrp[i]) {
          sigma33(j, k) += tmp * mat_x(i, j+1) * mat_x(i, k+1);
        }

        if (1 == mgrp[i]) {
          sigma23(j, k) += tmp2 * mat_x(i, j+1) * mat_x(i, k+1);
        }
      }
    }
  }

  // Sigma: take average
  for (j = 0; j < n_x + 1; j++) {
    for (k = j; k < n_x + 1; k++) {
      sigma11(j, k) /= n_pat;
      sigma11(k, j) =  sigma11(j, k);
    }
  }

  for (j = 0; j < n_x + 1; j++) {
    for (k = 0; k < n_x; k++) {
      sigma12(j, k) /= (n3[0] + n3[1]);
      sigma13(j, k) /= (n3[0] + n3[2]);
    }
  }

  for (j = 0; j < n_x; j++) {
    for (k = j; k < n_x; k++) {
      sigma22(j, k) /= (n3[0] + n3[1]);
      sigma33(j, k) /= (n3[0] + n3[2]);
      sigma23(j, k) /= n3[0];

      sigma22(k, j) = sigma22(j, k);
      sigma33(k, j) = sigma33(j, k);
      sigma23(k, j) = sigma23(j, k);
    }
  }

  // Sigma: put together
  for (j = 0; j < n_x + 1; j++) {
    for (k = 0; k < n_x + 1; k++) {
      sigma(j, k) = sigma11(j, k);
    }
  }

  for (j = 0; j < n_x; j++) {
    for (k = 0; k < n_x; k++) {
      sigma(n_x + 1 + j,     n_x + 1 + k    ) = sigma22(j, k);
      sigma(2 * n_x + 1 + j, 2 * n_x + 1 + k) = sigma33(j, k);
      sigma(n_x + 1 + j,     2 * n_x + 1 + k) = sigma23(j, k);
    }
  }

  // for (j = 0; j < n_x + 1; j++) {
  //   for (k = 0; k < n_x; k++) {
  //      sigma(j, n_x + 1 + k    ) = sigma12(j, k);
  //      sigma(j, 2 * n_x + 1 + k) = sigma13(j, k);
  //    }
  // }

  for (j = 0; j < 3 * n_x + 1; j++) {
    for (k = j; k < 3 * n_x + 1; k++) {
      sigma(k, j) = sigma(j, k);
    }
  }

  // return
  return(sigma);
}


// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
//                  TEMPORARY
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

//' Get Moment Constraints
//'
//' return N row k column matrix
//' Row:   subject
//' Column: Moment
//'
NumericMatrix c_ps_gmm_g_old (NumericVector beta,
                              NumericMatrix mat_grp_x,
                              bool att = false) {
  int n_pat  = mat_grp_x.nrow();
  int n_beta = beta.size();

  NumericMatrix rst(n_pat, 3 * n_beta - 2);
  NumericVector curx;

  int    i, k;
  double ex, z, r;
  double w;

  for (i = 0; i < n_pat; i++) {
    curx = mat_grp_x(i, _);

    if (curx(0) <= 1.0) {
      z = 1;
    } else {
      z = 0;
    }

    if (curx(0) <= 2.0) {
      r = 1;
    } else {
      r = 0;
    }

    ex = g_ex(beta, curx[Range(1, n_beta + 1)]);

    // set 1: grp 1 vs. grp 2:3
    if (att) {
      w =  z * ex * (1 - ex) - (1 - z) * ex * ex;
    } else {
      w =  z * (1 - ex) - (1 - z) * ex;
    }

    for (k = 0; k < n_beta; k++) {
      rst(i, k) = w * curx(k + 1);
    }

    // set 2: grp 2 vs. grp 3 wrt X
    // w =  (1 - z) * r - (1 - z) * (1 - r);
    // if (att) {
    //   w *=  ex / (1 - ex);
    // } else {
    //   w /=  (1 - ex);
    // }

    // set 2: grp 1 vs. grp 2 wrt X
    if (att) {
      w =  z - (1 - z) * r * ex / (1 - ex);
    } else {
      w =  z / ex - (1 - z) * r / (1 - ex);
    }

    for (k = 1; k < n_beta; k++) {
      rst(i, k + n_beta - 1) = w * curx(k + 1);
    }

    // set 3: grp 1 vs. grp 3 wrt X
    if (att) {
      w =  z - (1 - z) * (1 - r) * ex / (1 - ex);
    } else {
      w =  z / ex - (1 - z) * (1 - r) / (1 - ex);
    }

    for (k = 1; k < n_beta; k++) {
      rst(i, k + 2 * n_beta - 2) = w * curx(k + 1);
    }

  }

  // return
  return(rst);
}

//' Get Derivative of Moment Constraints
//'
//' return k row p column matrix
//' Row:    1..k moment functions
//' Column: 1..p coefficients
//'
//'
// [[Rcpp::export]]
NumericMatrix c_ps_gmm_dg (NumericVector beta,
                          NumericMatrix mat_grp_x,
                          bool att = false) {
  int n_pat  = mat_grp_x.nrow();
  int n_beta = beta.size();
  NumericMatrix rst(2 * n_beta - 1, n_beta);
  NumericVector curx;

  int i, k, l;
  double ex, z, r;
  double w;

  for (i = 0; i < n_pat; i++) {
    curx = mat_grp_x(i, _);

    if (curx(0) <= 1.0) {
      z = 1;
    } else {
      z = 0;
    }

    if (curx(0) <= 2.0) {
      r = 1;
    } else {
      r = 0;
    }

    ex = g_ex(beta, curx[Range(1, n_beta + 1)]);

    // set 1: grp 1 vs. grp 2:3
    if (att) {
      w = z * (1 - 2 * ex) - 2 * (1 - z) * ex;
    } else {
      w = - z - (1 - z);
    }

    for (k = 0; k < n_beta; k++) {
      for (l = 0; l < n_beta; l++) {
        rst(k, l) += w * ex * (1 - ex) * curx(k + 1) * curx(l + 1);
      }
    }

    // set 2: grp 2 vs. grp 3 wrt X
    w =  (1 - z) * r - (1 - z) * (1 - r);
    // w * exp(xb)
    w *= ex / (1 - ex);

    for (k = 1; k < n_beta; k++) {
      for (l = 0; l < n_beta; l++) {
        rst(k + n_beta - 1, l) += w * curx(k + 1) * curx(l + 1);
      }
    }

    // set 2: grp 1 vs. grp 2
    // if (att) {
    //   w = z * (1 - 2 * ex) - 2 * (1 - z) * r * ex;
    // } else {
    //   w = - z - (1 - z) * r;
    // }

    // for (k = 0; k < n_beta; k++) {
    //   for (l = 0; l < n_beta; l++) {
    //     rst(k + n_beta, l) += w * ex * (1 - ex) * curx(k + 1) * curx(l + 1);
    //   }
    // }

    // set 3: grp 1 vs. grp 3
    // if (att) {
    //   w = z * (1 - 2 * ex) - 2 * (1 - z) * (1 - r) * ex;
    // } else {
    //   w = - z - (1 - z) * (1 - r);
    // }

    // for (k = 0; k < n_beta; k++) {
    //   for (l = 0; l < n_beta; l++) {
    //     rst(k + 2 * n_beta, l) +=
    //       w * ex * (1 - ex) * curx(k + 1) * curx(l + 1);
    //   }
    // }
  }

  // take average
  for (k = 0; k < rst.nrow(); k++) {
    for (l = 0; l < rst.ncol(); l++) {
      rst(k, l) /= n_pat;
    }
  }

  // return
  return(rst);
}
