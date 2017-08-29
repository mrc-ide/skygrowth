//~ #include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;



vec mh_sample_ne1_2( vec fwdnevec, vec fwdnco , double tau , vec prop_sigma, double dh, mat lterms);
RcppExport SEXP phyloma_mh_sample_ne1_2( SEXP fwdnevecSEXP, SEXP fwdncoSEXP,  SEXP tauSEXP, SEXP prop_sigmaSEXP, SEXP dhSEXP , SEXP ltermsSEXP){
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type nevec(fwdnevecSEXP);
    Rcpp::traits::input_parameter< vec >::type nco(fwdncoSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< vec >::type prop_sigma(prop_sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type dh(dhSEXP);
    Rcpp::traits::input_parameter< mat >::type lterms(ltermsSEXP);
    rcpp_result_gen = Rcpp::wrap(mh_sample_ne1_2(nevec, nco, tau, prop_sigma, dh, lterms));
    return rcpp_result_gen;
END_RCPP
}

vec mh_sample_ne2_2( vec fwdnevec, vec fwdnco,  double tau , vec prop_sigma, double dh, vec zxb, mat lterms);
RcppExport SEXP phyloma_mh_sample_ne2_2( SEXP fwdnevecSEXP, SEXP fwdncoSEXP, SEXP tauSEXP, SEXP prop_sigmaSEXP, SEXP dhSEXP, SEXP zxbSEXP, SEXP ltermsSEXP ){
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type nevec(fwdnevecSEXP);
    Rcpp::traits::input_parameter< vec >::type nco(fwdncoSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< vec >::type prop_sigma(prop_sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type dh(dhSEXP);
    Rcpp::traits::input_parameter< vec >::type zxb(zxbSEXP);
    Rcpp::traits::input_parameter< mat >::type lterms(ltermsSEXP);
    rcpp_result_gen = Rcpp::wrap(mh_sample_ne2_2(nevec, nco,  tau, prop_sigma, dh, zxb, lterms));
    return rcpp_result_gen;
END_RCPP
}

