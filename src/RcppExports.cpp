#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// insert_sorted_cpp
NumericVector insert_sorted_cpp(NumericVector sorted_x, double value);
RcppExport SEXP _RepsampleR_insert_sorted_cpp(SEXP sorted_xSEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type sorted_x(sorted_xSEXP);
    Rcpp::traits::input_parameter< double >::type value(valueSEXP);
    rcpp_result_gen = Rcpp::wrap(insert_sorted_cpp(sorted_x, value));
    return rcpp_result_gen;
END_RCPP
}
// ks2_d_with_candidate_cpp
double ks2_d_with_candidate_cpp(NumericVector sample_sorted, double candidate, NumericVector pop_sorted);
RcppExport SEXP _RepsampleR_ks2_d_with_candidate_cpp(SEXP sample_sortedSEXP, SEXP candidateSEXP, SEXP pop_sortedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type sample_sorted(sample_sortedSEXP);
    Rcpp::traits::input_parameter< double >::type candidate(candidateSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pop_sorted(pop_sortedSEXP);
    rcpp_result_gen = Rcpp::wrap(ks2_d_with_candidate_cpp(sample_sorted, candidate, pop_sorted));
    return rcpp_result_gen;
END_RCPP
}
// ks1_norm_d_with_candidate_cpp
double ks1_norm_d_with_candidate_cpp(NumericVector sample_sorted, double candidate, double mean, double sd);
RcppExport SEXP _RepsampleR_ks1_norm_d_with_candidate_cpp(SEXP sample_sortedSEXP, SEXP candidateSEXP, SEXP meanSEXP, SEXP sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type sample_sorted(sample_sortedSEXP);
    Rcpp::traits::input_parameter< double >::type candidate(candidateSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    rcpp_result_gen = Rcpp::wrap(ks1_norm_d_with_candidate_cpp(sample_sorted, candidate, mean, sd));
    return rcpp_result_gen;
END_RCPP
}
// ks1_cdf_d_with_candidate_cpp
double ks1_cdf_d_with_candidate_cpp(NumericVector sample_sorted, NumericVector sample_cdf_sorted, double candidate, double candidate_cdf);
RcppExport SEXP _RepsampleR_ks1_cdf_d_with_candidate_cpp(SEXP sample_sortedSEXP, SEXP sample_cdf_sortedSEXP, SEXP candidateSEXP, SEXP candidate_cdfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type sample_sorted(sample_sortedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sample_cdf_sorted(sample_cdf_sortedSEXP);
    Rcpp::traits::input_parameter< double >::type candidate(candidateSEXP);
    Rcpp::traits::input_parameter< double >::type candidate_cdf(candidate_cdfSEXP);
    rcpp_result_gen = Rcpp::wrap(ks1_cdf_d_with_candidate_cpp(sample_sorted, sample_cdf_sorted, candidate, candidate_cdf));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RepsampleR_insert_sorted_cpp", (DL_FUNC) &_RepsampleR_insert_sorted_cpp, 2},
    {"_RepsampleR_ks2_d_with_candidate_cpp", (DL_FUNC) &_RepsampleR_ks2_d_with_candidate_cpp, 3},
    {"_RepsampleR_ks1_norm_d_with_candidate_cpp", (DL_FUNC) &_RepsampleR_ks1_norm_d_with_candidate_cpp, 4},
    {"_RepsampleR_ks1_cdf_d_with_candidate_cpp", (DL_FUNC) &_RepsampleR_ks1_cdf_d_with_candidate_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_RepsampleR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
