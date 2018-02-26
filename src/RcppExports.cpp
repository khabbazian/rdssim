// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// adj2list
Rcpp::List adj2list(SEXP X_, std::string acType);
RcppExport SEXP rdssim_adj2list(SEXP X_SEXP, SEXP acTypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< std::string >::type acType(acTypeSEXP);
    rcpp_result_gen = Rcpp::wrap(adj2list(X_, acType));
    return rcpp_result_gen;
END_RCPP
}

// rdssim_cpp
Rcpp::NumericMatrix rdssim_cpp(Rcpp::List rcpp_adjList, Rcpp::List rcpp_acAdjList, std::string referralType, bool wReplacement, int nSamples, int nReferrals, int seedNode, int rseed);
RcppExport SEXP rdssim_rdssim_cpp(SEXP rcpp_adjListSEXP, SEXP rcpp_acAdjListSEXP, SEXP referralTypeSEXP, SEXP wReplacementSEXP, SEXP nSamplesSEXP, SEXP nReferralsSEXP, SEXP seedNodeSEXP, SEXP rseedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type rcpp_adjList(rcpp_adjListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type rcpp_acAdjList(rcpp_acAdjListSEXP);
    Rcpp::traits::input_parameter< std::string >::type referralType(referralTypeSEXP);
    Rcpp::traits::input_parameter< bool >::type wReplacement(wReplacementSEXP);
    Rcpp::traits::input_parameter< int >::type nSamples(nSamplesSEXP);
    Rcpp::traits::input_parameter< int >::type nReferrals(nReferralsSEXP);
    Rcpp::traits::input_parameter< int >::type seedNode(seedNodeSEXP);
    Rcpp::traits::input_parameter< int >::type rseed(rseedSEXP);
    rcpp_result_gen = Rcpp::wrap(rdssim_cpp(rcpp_adjList, rcpp_acAdjList, referralType, wReplacement, nSamples, nReferrals, seedNode, rseed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"rdssim_adj2list", (DL_FUNC) &rdssim_adj2list, 2},
    {"rdssim_rdssim_cpp", (DL_FUNC) &rdssim_rdssim_cpp, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_rdssim(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
