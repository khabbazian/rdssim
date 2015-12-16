#include <RcppEigen.h>
#include <Rcpp.h>
#include "types.h"

//using namespace Rcpp;

// adj2list
AdjList adj2list(SEXP X_);
RcppExport SEXP rdssim_adj2list(SEXP X_SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type X_(X_SEXP);
    __result = Rcpp::wrap(adj2list(X_));
    return __result;
END_RCPP
}
// rdssimMarkov
Matrix rdssimMarkov(Rcpp::List rcpp_adjlist, string rType, int nSamples, int nReferrals, int seedNode, int rseed);
RcppExport SEXP rdssim_rdssimMarkov(SEXP rcpp_adjlistSEXP, SEXP rTypeSEXP, SEXP nSamplesSEXP, SEXP nReferralsSEXP, SEXP seedNodeSEXP, SEXP rseedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::List >::type rcpp_adjlist(rcpp_adjlistSEXP);
    Rcpp::traits::input_parameter< string >::type rType(rTypeSEXP);
    Rcpp::traits::input_parameter< int >::type nSamples(nSamplesSEXP);
    Rcpp::traits::input_parameter< int >::type nReferrals(nReferralsSEXP);
    Rcpp::traits::input_parameter< int >::type seedNode(seedNodeSEXP);
    Rcpp::traits::input_parameter< int >::type rseed(rseedSEXP);
    __result = Rcpp::wrap(rdssimMarkov(rcpp_adjlist, rType, nSamples, nReferrals, seedNode, rseed));
    return __result;
END_RCPP
}
