#include <armadillo>
#include <RcppArmadillo.h>
#include <boost/dynamic_bitset.hpp>
#include <vector>

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
std::vector< std::string > signature_mapper(arma::mat X) {
    string s;
    vector<string> signatures;
    
    // if the min is not 0, then find the min and round it in case there
    // are problems with presicion;
    // this should be changed if the cpm normalisation will become standard in
    // scater package - at the moment cpm convert all zeros to non-zeros...
    float zero = ceilf(X.min() * 100) / 100;
    for(int j = 0; j < X.n_cols; j++) {
        // construct a signature
        s = "";
        for(int i = 0; i < X.n_rows; i++) {
            if(X(i, j) > zero) {
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
std::vector<std::string> sort_lexic(std::vector<std::string> v) {
    std::sort(v.begin(), v.end());
    return(v);
}

// [[Rcpp::export]]
arma::vec signatures_order(std::vector< std::string > signatures) {
    int i, j;
    int k;
    arma::vec res = zeros<vec>(signatures.size());
    std::string s;
    for(i = 0; i < signatures.size(); i++) {
        s = signatures[i];
        k = 0;
        for(j = 0; j < s.length(); j++) {
            if(s[j] == '1') {
                k++;
            }
        }
        res(i) = k;
    }
    return(res);
}

// [[Rcpp::export]]
std::vector< std::string > signature_mapper_fjlt(arma::mat X, arma::vec means) {
    string s;
    vector<string> signatures;
    
    for(int j = 0; j < X.n_cols; j++) {
        // construct a signature
        s = "";
        for(int i = 0; i < X.n_rows; i++) {
            if(X(i, j) >= means(i)) {
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
        buckets(i) = b;
        for(j = i + 1; j < signatures.size(); j++) {
            if(buckets(j) == 0) {
                // create a bitset from string j
                boost::dynamic_bitset<> b2(signatures[j]);
                // XOR the bitsets, when bits are equal, the result is 0, otherwise 1
                b3 = b1 ^ b2;
                // count all 1 in the resulting bitset and subtract from the size
                // this gives the total number of matching bits between the bitsets
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

