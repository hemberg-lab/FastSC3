#include <armadillo>
#include <RcppArmadillo.h>
#include <boost/dynamic_bitset.hpp>

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec get_span(arma::mat X) {
    int i;
    arma::rowvec v;
    arma::vec span(X.n_rows);
    for(i = 0; i < X.n_rows; i++) {
        v = X.row(i);
        span(i) = max(v) - min(v);
    }
    return(span);
}

// [[Rcpp::export]]
arma::vec get_hyperplane(arma::vec span, int M) {
    span = cumsum(span/sum(span));
    arma::vec rnd = randu<vec>(M);
    arma::vec inds(M);
    arma::vec span1;
    for(int i = 0; i < M; i++) {
        span1 = span.elem( find( span < rnd(i) ) );
        inds(i) = span1.size() + 1;
    }
    return(inds);
}

// [[Rcpp::export]]
arma::vec get_thresholds(arma::mat X, int bin_num) {
    arma::vec bin(bin_num);
    arma::rowvec v;
    arma::vec v1;
    arma::vec span(X.n_rows);
    arma::vec threshold(X.n_rows);
    for(int i = 0; i < X.n_rows; i++) {
        v = X.row(i);
        span(i) = max(v) - min(v);
        for(int j = 0; j < bin_num; j++) {
            v1 = v.elem( find(v >= (min(v) + j * span(i) / bin_num) && v <= (min(v) + j * span(i) / bin_num)) );
            bin(j) = v1.size();
        }
        threshold(i) = min(v) + bin.index_min() * span(i) / bin_num;
    }
    return(threshold);
}

// [[Rcpp::export]]
std::vector< std::string > signature_mapper(arma::mat X, int M, int bin_num) {
    arma::vec hyperplane;
    string s;
    vector<string> signatures;
    // calculate all spans
    arma::vec span = get_span(X);
    // calculate all thresholds
    arma::vec thresholds = get_thresholds(X, bin_num);
    for(int j = 0; j < X.n_cols; j++) {
        // define a hyperplane
        hyperplane = get_hyperplane(span, M);
        // construct a signature
        s = "";
        for(int i = 0; i < M; i++) {
            if(X(hyperplane(i), j) <= thresholds(hyperplane(i))) {
                s += '1';
            } else {
                s += '0';
            }
        }
        signatures.push_back(s);
    }
    return(signatures);
}

// [[Rcpp::export]]
arma::vec get_buckets(std::vector< std::string > signatures, int P) {
    int i, j, b = 1;
    boost::dynamic_bitset<> b3;
    arma::vec buckets = zeros<vec>(signatures.size());
    while(min(buckets) == 0) {
        // find first unassigned column
        i = 0;
        while(buckets(i) != 0) {
            i++;
        }
        boost::dynamic_bitset<> b1(signatures[i]);
        for(j = i; j < signatures.size(); j++) {
            if(buckets(j) == 0) {
                // create a bitset from string j
                boost::dynamic_bitset<> b2(signatures[j]);
                // XOR the bitsets, when bits are equal, the result is 0, otherwise 1
                b3 = b1 ^ b2;
                // count all 1 in the resulting bitset and subtract from the size
                // this give the total number of matching bits between the bitsets
                if( ( b3.size() - b3.count() ) >= P ) {
                    buckets(j) = b;
                }
            }
        }
        // move to the next bucket
        b++;
    }
    return(buckets);
}

