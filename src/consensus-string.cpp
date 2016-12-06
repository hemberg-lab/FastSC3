#include <Rcpp.h>
#include <boost/dynamic_bitset.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
std::string get_consensus_string(std::vector< std::string > signatures, double threshold) {
    int i, j;
    
    std::string res = "";
    
    int zeros = 0, ones = 0;
    
    for(i = 0; i < signatures[0].size(); i++) {
        for(j = 0; j < signatures.size(); j++) {
            if(signatures[j].at(i) == '0') {
                zeros++;
            } else {
                ones++;
            }
        }
        if(ones > threshold * (zeros + ones)) {
            res += "1";
        } else if(zeros > threshold * (zeros + ones)) {
            res += "0";
        } else {
            res += "_";
        }
        zeros = 0;
        ones = 0;
    }
    return(res);
}

// [[Rcpp::export]]
std::vector< int > compare_signatures(std::string bsig, std::vector< std::string > signatures) {
    boost::dynamic_bitset<> bs(bsig);
    boost::dynamic_bitset<> sxor;
    std::vector< int > res;
    
    for(int i = 0; i < signatures.size(); i++) {
        boost::dynamic_bitset<> cs(signatures[i]);
        sxor = bs ^ cs;
        res.push_back(sxor.size() - sxor.count());
    }
    
    return(res);
}
