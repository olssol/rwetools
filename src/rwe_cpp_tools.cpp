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

  int i, j, k, start;

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


//' Get Moment Constraints
//'
//' return N row k column matrix
//' Row:   subject
//' Column: Moment
//'
//'
// [[Rcpp::export]]
NumericMatrix c_ps_gmm_g (NumericVector beta,
                          NumericMatrix mat_grp_x,
                          bool att = false) {
  int n_pat  = mat_grp_x.nrow();
  int n_beta = beta.size();

  NumericMatrix rst(n_pat, 3 * n_beta - 2);
  NumericVector curx;

  int    i, j, k;
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

  int i, j, k, l;
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
