// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rcpp_distance
NumericVector rcpp_distance(int iJob, int nJobs, List sparseMatrix, const char* method, bool weighted);
RcppExport SEXP rbiom_rcpp_distance(SEXP iJobSEXP, SEXP nJobsSEXP, SEXP sparseMatrixSEXP, SEXP methodSEXP, SEXP weightedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type iJob(iJobSEXP);
    Rcpp::traits::input_parameter< int >::type nJobs(nJobsSEXP);
    Rcpp::traits::input_parameter< List >::type sparseMatrix(sparseMatrixSEXP);
    Rcpp::traits::input_parameter< const char* >::type method(methodSEXP);
    Rcpp::traits::input_parameter< bool >::type weighted(weightedSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_distance(iJob, nJobs, sparseMatrix, method, weighted));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_unifrac
NumericVector rcpp_unifrac(int iJob, int nJobs, List sparseMatrix, List tree, bool weighted);
RcppExport SEXP rbiom_rcpp_unifrac(SEXP iJobSEXP, SEXP nJobsSEXP, SEXP sparseMatrixSEXP, SEXP treeSEXP, SEXP weightedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type iJob(iJobSEXP);
    Rcpp::traits::input_parameter< int >::type nJobs(nJobsSEXP);
    Rcpp::traits::input_parameter< List >::type sparseMatrix(sparseMatrixSEXP);
    Rcpp::traits::input_parameter< List >::type tree(treeSEXP);
    Rcpp::traits::input_parameter< bool >::type weighted(weightedSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_unifrac(iJob, nJobs, sparseMatrix, tree, weighted));
    return rcpp_result_gen;
END_RCPP
}
