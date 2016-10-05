#include <armadillo>
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

//' Calculate approximate eigenvectors
//' 
//' Nystrom spectral shifting method
//' 
//' @param K input SPSD matrix
//' @param p number of selected columns as fraction of ncol(K)
//' 
//' @return approximate eigenvectors
//' 
// [[Rcpp::export]]
arma::mat ssNystrom(arma::mat K, int r) {
    // create a vector v = 1:N
    std::vector<int> v;
    for(int i = 0; i < K.n_cols; i++) {
        v.push_back(i);
    }
    
    // permute v and write to inds
    arma::vec N_vec = conv_to<vec>::from(v);
    arma::vec inds = shuffle( N_vec );
    
    // subset inds
    inds = inds.subvec(0, r - 1);
    
    // subset K using inds
    NumericVector idx = NumericVector(inds.begin(), inds.end());
    arma::uvec idx1 = as<uvec>(idx);
    arma::mat C = K.cols(idx1);
    
    // compute Moore–Penrose pseudoinverse
    arma::mat Ci = pinv( C );
    
    // compute U (rxr)
    arma::mat U = Ci * K * Ci.t();
    
    // SVD of C
    mat Uc;
    vec Ec;
    mat Vc;
    svd_econ(Uc,Ec,Vc,C, "both", "std");
    
    // Calculate S
    mat S = diagmat(Ec) * Vc * U * Vc.t() * diagmat(Ec).t();

    // Eigenvalue decomposition of S
    cx_vec Ls;
    cx_mat Us;

    eig_gen(Ls, Us, S);

    // compute approximate eigenvectors
    mat Us_conv = conv_to<mat>::from(Us);
    arma::mat apprEigen = Uc * Us_conv;
    
    return(apprEigen);
}

