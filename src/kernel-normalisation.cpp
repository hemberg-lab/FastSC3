#include <armadillo>
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat calculate_P1(arma::mat X) {
    int n = X.n_cols;
    arma::mat P1 = X + ( (1/n + cumsum(X) / pow(n, 2)) * eye(n, n) - 1/n * X ) - 1/n * ones(n, n) * X;
    return(P1);
}

// [[Rcpp::export]]
arma::mat normalise_kernel(arma::mat K) {
    arma::mat X = calculate_P1(K);
    while(X.min() < 0.0) {
        X.elem( find(X < 0.0) ).zeros();
        X = calculate_P1(X);
    }
    return(X);
}
