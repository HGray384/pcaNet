// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// bpcaNet
List bpcaNet(arma::mat myMat, int N, int D, arma::uvec hidden, arma::uvec numberOfNonNAvaluesInEachCol, arma::uvec nomissIndex, arma::uvec missIndex, int nMissing, int nPcs, double threshold, int maxIterations);
RcppExport SEXP _ppcaNet_bpcaNet(SEXP myMatSEXP, SEXP NSEXP, SEXP DSEXP, SEXP hiddenSEXP, SEXP numberOfNonNAvaluesInEachColSEXP, SEXP nomissIndexSEXP, SEXP missIndexSEXP, SEXP nMissingSEXP, SEXP nPcsSEXP, SEXP thresholdSEXP, SEXP maxIterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type myMat(myMatSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type hidden(hiddenSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type numberOfNonNAvaluesInEachCol(numberOfNonNAvaluesInEachColSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type nomissIndex(nomissIndexSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type missIndex(missIndexSEXP);
    Rcpp::traits::input_parameter< int >::type nMissing(nMissingSEXP);
    Rcpp::traits::input_parameter< int >::type nPcs(nPcsSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< int >::type maxIterations(maxIterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(bpcaNet(myMat, N, D, hidden, numberOfNonNAvaluesInEachCol, nomissIndex, missIndex, nMissing, nPcs, threshold, maxIterations));
    return rcpp_result_gen;
END_RCPP
}
// ppcaNet
List ppcaNet(arma::mat myMat, int N, int D, arma::mat W, arma::uvec hidden, int nMissing, int nPcs, double threshold, int maxIterations);
RcppExport SEXP _ppcaNet_ppcaNet(SEXP myMatSEXP, SEXP NSEXP, SEXP DSEXP, SEXP WSEXP, SEXP hiddenSEXP, SEXP nMissingSEXP, SEXP nPcsSEXP, SEXP thresholdSEXP, SEXP maxIterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type myMat(myMatSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type hidden(hiddenSEXP);
    Rcpp::traits::input_parameter< int >::type nMissing(nMissingSEXP);
    Rcpp::traits::input_parameter< int >::type nPcs(nPcsSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< int >::type maxIterations(maxIterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(ppcaNet(myMat, N, D, W, hidden, nMissing, nPcs, threshold, maxIterations));
    return rcpp_result_gen;
END_RCPP
}
// ppcaSensible
List ppcaSensible(const arma::mat Y, arma::mat W, double v, const double traceS, const int MaxIter, const double TolFun, const double TolX);
RcppExport SEXP _ppcaNet_ppcaSensible(SEXP YSEXP, SEXP WSEXP, SEXP vSEXP, SEXP traceSSEXP, SEXP MaxIterSEXP, SEXP TolFunSEXP, SEXP TolXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< const double >::type traceS(traceSSEXP);
    Rcpp::traits::input_parameter< const int >::type MaxIter(MaxIterSEXP);
    Rcpp::traits::input_parameter< const double >::type TolFun(TolFunSEXP);
    Rcpp::traits::input_parameter< const double >::type TolX(TolXSEXP);
    rcpp_result_gen = Rcpp::wrap(ppcaSensible(Y, W, v, traceS, MaxIter, TolFun, TolX));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP ppcaNet_bpcaNet(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP ppcaNet_ppcaNet(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_ppcaNet_bpcaNet", (DL_FUNC) &_ppcaNet_bpcaNet, 11},
    {"_ppcaNet_ppcaNet", (DL_FUNC) &_ppcaNet_ppcaNet, 9},
    {"_ppcaNet_ppcaSensible", (DL_FUNC) &_ppcaNet_ppcaSensible, 7},
    {"ppcaNet_bpcaNet", (DL_FUNC) &ppcaNet_bpcaNet, 11},
    {"ppcaNet_ppcaNet", (DL_FUNC) &ppcaNet_ppcaNet,  9},
    {NULL, NULL, 0}
};

RcppExport void R_init_ppcaNet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
