// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// constructP
arma::mat constructP(int p, int k, int d, int n);
RcppExport SEXP FastSC3_constructP(SEXP pSEXP, SEXP kSEXP, SEXP dSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    __result = Rcpp::wrap(constructP(p, k, d, n));
    return __result;
END_RCPP
}
// constructH
arma::mat constructH(int d);
RcppExport SEXP FastSC3_constructH(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    __result = Rcpp::wrap(constructH(d));
    return __result;
END_RCPP
}
// constructD
arma::rowvec constructD(int d);
RcppExport SEXP FastSC3_constructD(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    __result = Rcpp::wrap(constructD(d));
    return __result;
END_RCPP
}
// calculateFJLT
arma::mat calculateFJLT(arma::mat x, int p, int k, int d, int n);
RcppExport SEXP FastSC3_calculateFJLT(SEXP xSEXP, SEXP pSEXP, SEXP kSEXP, SEXP dSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    __result = Rcpp::wrap(calculateFJLT(x, p, k, d, n));
    return __result;
END_RCPP
}
