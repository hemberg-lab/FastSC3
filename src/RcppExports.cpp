// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// pad_matrix_with_zeros
arma::mat pad_matrix_with_zeros(arma::mat x, int d);
RcppExport SEXP FastSC3_pad_matrix_with_zeros(SEXP xSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(pad_matrix_with_zeros(x, d));
    return rcpp_result_gen;
END_RCPP
}
// constr_P
arma::mat constr_P(int p, int k, int d, int n);
RcppExport SEXP FastSC3_constr_P(SEXP pSEXP, SEXP kSEXP, SEXP dSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(constr_P(p, k, d, n));
    return rcpp_result_gen;
END_RCPP
}
// constr_H
arma::mat constr_H(int d);
RcppExport SEXP FastSC3_constr_H(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(constr_H(d));
    return rcpp_result_gen;
END_RCPP
}
// constr_D
arma::rowvec constr_D(int d);
RcppExport SEXP FastSC3_constr_D(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(constr_D(d));
    return rcpp_result_gen;
END_RCPP
}
// calc_fjlt
arma::mat calc_fjlt(arma::mat x, int p, int k, int d, int n);
RcppExport SEXP FastSC3_calc_fjlt(SEXP xSEXP, SEXP pSEXP, SEXP kSEXP, SEXP dSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_fjlt(x, p, k, d, n));
    return rcpp_result_gen;
END_RCPP
}
// calculate_P1
arma::mat calculate_P1(arma::mat X);
RcppExport SEXP FastSC3_calculate_P1(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_P1(X));
    return rcpp_result_gen;
END_RCPP
}
// normalise_kernel
arma::mat normalise_kernel(arma::mat K);
RcppExport SEXP FastSC3_normalise_kernel(SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(normalise_kernel(K));
    return rcpp_result_gen;
END_RCPP
}
// signature_mapper
std::vector< std::string > signature_mapper(arma::mat X);
RcppExport SEXP FastSC3_signature_mapper(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(signature_mapper(X));
    return rcpp_result_gen;
END_RCPP
}
// get_buckets
arma::vec get_buckets(std::vector< std::string > signatures, int P);
RcppExport SEXP FastSC3_get_buckets(SEXP signaturesSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector< std::string > >::type signatures(signaturesSEXP);
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(get_buckets(signatures, P));
    return rcpp_result_gen;
END_RCPP
}
// calc_delta
double calc_delta(arma::mat K, int k);
RcppExport SEXP FastSC3_calc_delta(SEXP KSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_delta(K, k));
    return rcpp_result_gen;
END_RCPP
}
// ssNystrom
arma::mat ssNystrom(arma::mat K, int c);
RcppExport SEXP FastSC3_ssNystrom(SEXP KSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(ssNystrom(K, c));
    return rcpp_result_gen;
END_RCPP
}
