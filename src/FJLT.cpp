#include <armadillo>
#include <RcppArmadillo.h>
#include <random>

using namespace std;
using namespace Rcpp;
using namespace arma;

//' Construct the P matrix for the FJLT transform
//' 
//' @param epsilon error tolerance parameter
//' @param p the norm
//' @param k dimension we reduce to
//' @param d dimension of the input matrix
//' @param n number of cells (columns) in the input matrix
//' 
//' @return matrix P of the FJLT transform
//' 
// [[Rcpp::export]]
arma::mat constr_P(int p, int k, int d, int n) {
    unsigned int i, j;
    
    // define q
    // p is 2, therefore we can remove the first term from the criteria
    double term = pow(log(n), p)/d;
    double q = min(term, 1.0);
    
    // initialize P
    rowvec tmp = zeros<rowvec>(3);
    mat P = zeros<mat>(k, d);
    
    // assign values to P
    for (i = 0; i < k; i++) {
        for (j = 0; j < d; j++) {
            if(randu() <= q) {
                P(i, j) = randn() * sqrt(1/q);
            }
        }
    }
    return(P);
}

//' Construct the H (Hadamard) matrix for the FJLT transform
//' 
//' @param d dimension of the input matrix
//' 
//' @return matrix H of the FJLT transform
//' 
// [[Rcpp::export]]
arma::mat constr_H(int d) {
    // initialize H
    mat H = ones<mat>(1, 1);
    while(d > 0) {
        H = join_rows(join_cols(H, H), join_cols(H, H*(-1)));
        d--;
    }
    return(H);
}

//' Construct the D matrix for the FJLT transform
//' 
//' @param d dimension of the input matrix
//' 
//' @return matrix D of the FJLT transform
//' 
// [[Rcpp::export]]
arma::rowvec constr_D(int d) {
    int i;
    rowvec D = zeros<rowvec>(d);
    for (i = 0; i < d; i++) {
        if(randu() <= 0.5) {
            D(i) = 1;
        } else {
            D(i) = -1;
        }
    }
    return(D);
}

//' Calculate the FJLT transform
//' 
//' Fast-Johnson-Lindenstrauss-Transform
//' 
//' @param epsilon error tolerance parameter
//' @param p the norm
//' @param k dimension we reduce to
//' @param d dimension of the input matrix
//' @param n number of cells (columns) in the input matrix
//' 
//' @return FJLT transform
//' 
// [[Rcpp::export]]
arma::mat calc_fjlt(arma::mat x, int p, int k, int d, int n) {

    mat P = constr_P(p, k, d, n);
    mat H = constr_H(log2(d))/sqrt(d);
    rowvec D = constr_D(d);
    
    // D is a vector of diagonal elements therefore we can just
    // multiply each column of H by the elements of D
    H.each_row() %= D;

    mat res = P * H * x;

    return(res);
}

