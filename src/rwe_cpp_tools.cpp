#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::init]]
void my_package_init(DllInfo *dll) {
  // initialization code here
  R_useDynamicSymbols(dll, TRUE);
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
