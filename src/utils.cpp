#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::mat t_cpp(arma::mat X) {
    return(X.t());
}
