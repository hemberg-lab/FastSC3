#include <armadillo>
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

//' Calculate the FJLT transform
//' 
//' Fast-Johnson-Lindenstrauss-Transform
//' 
//' @param x input matrix
//' @param d dimension to pad to
//' 
//' @return FJLT transform
//' 

// [[Rcpp::export]]
arma::mat CalculateApproximateLaplacianEigenvectors(arma::mat L, double p, int d, double delta = 0.9985) {
    int N = L.n_cols;
    int r = round(N * p);
    
    // create a vector v = 1:N
    std::vector<int> v;
    for(int i = 0; i < N; i++) {
        v.push_back(i);
    }
    
    // permute v and write to sampleInds
    arma::vec N_vec = conv_to<vec>::from(v);
    arma::vec sampleInds = shuffle( N_vec );
    
    // subset sampleInds
    sampleInds = sampleInds.subvec(0, r - 1);
    
    // subset L using sampleInds
    NumericVector idx = NumericVector(sampleInds.begin(), sampleInds.end());
    arma::uvec idx1 = as<uvec>(idx);
    arma::mat C = L.cols(idx1);
    
    // compute Moore–Penrose pseudoinverse
    arma::mat Ci = pinv( C );
    
    // compute U
    arma::mat U = Ci * L * Ci.t();
    
    // compute modified SS-Nystrom
    arma::mat K = C * U * C.t();
    
    return(K);
}

// // [[Rcpp::export]]
// arma::vec test(int N) {
//     
//     return(sampleInds.subvec(0, 5));
// }

//' Calculate the FJLT transform
//' 
//' Fast-Johnson-Lindenstrauss-Transform
//' 
//' @param x input matrix
//' @param d dimension to pad to
//' 
//' @return FJLT transform
//' 
// [[Rcpp::export]]
arma::mat AdaptiveSampling(arma::mat x, int d) {
    x.insert_rows( x.n_rows, d - x.n_rows );
    return(x);
}
