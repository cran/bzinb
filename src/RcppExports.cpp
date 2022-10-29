// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// em
List em(NumericVector& param2, IntegerVector& xvec, IntegerVector& yvec, IntegerVector& freq, int& n, int& se, int& maxiter, double& tol, int showFlag, int bnb);
RcppExport SEXP _bzinb_em(SEXP param2SEXP, SEXP xvecSEXP, SEXP yvecSEXP, SEXP freqSEXP, SEXP nSEXP, SEXP seSEXP, SEXP maxiterSEXP, SEXP tolSEXP, SEXP showFlagSEXP, SEXP bnbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type param2(param2SEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type xvec(xvecSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type yvec(yvecSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type freq(freqSEXP);
    Rcpp::traits::input_parameter< int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< int& >::type se(seSEXP);
    Rcpp::traits::input_parameter< int& >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type showFlag(showFlagSEXP);
    Rcpp::traits::input_parameter< int >::type bnb(bnbSEXP);
    rcpp_result_gen = Rcpp::wrap(em(param2, xvec, yvec, freq, n, se, maxiter, tol, showFlag, bnb));
    return rcpp_result_gen;
END_RCPP
}
// dBvZINB_Expt
void dBvZINB_Expt(int& x, int& y, int& freq, double& a0, double& a1, double& a2, double& b1, double& b2, double& p1, double& p2, double& p3, double& p4, NumericVector& expt, NumericVector& s_i, NumericVector& info, int se, int bnb);
RcppExport SEXP _bzinb_dBvZINB_Expt(SEXP xSEXP, SEXP ySEXP, SEXP freqSEXP, SEXP a0SEXP, SEXP a1SEXP, SEXP a2SEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP p1SEXP, SEXP p2SEXP, SEXP p3SEXP, SEXP p4SEXP, SEXP exptSEXP, SEXP s_iSEXP, SEXP infoSEXP, SEXP seSEXP, SEXP bnbSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int& >::type y(ySEXP);
    Rcpp::traits::input_parameter< int& >::type freq(freqSEXP);
    Rcpp::traits::input_parameter< double& >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< double& >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< double& >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< double& >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< double& >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< double& >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< double& >::type p2(p2SEXP);
    Rcpp::traits::input_parameter< double& >::type p3(p3SEXP);
    Rcpp::traits::input_parameter< double& >::type p4(p4SEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type expt(exptSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type s_i(s_iSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type info(infoSEXP);
    Rcpp::traits::input_parameter< int >::type se(seSEXP);
    Rcpp::traits::input_parameter< int >::type bnb(bnbSEXP);
    dBvZINB_Expt(x, y, freq, a0, a1, a2, b1, b2, p1, p2, p3, p4, expt, s_i, info, se, bnb);
    return R_NilValue;
END_RCPP
}
// dBvZINB_Expt_vec
void dBvZINB_Expt_vec(IntegerVector& xvec, IntegerVector& yvec, IntegerVector& freq, int& n, double& a0, double& a1, double& a2, double& b1, double& b2, double& p1, double& p2, double& p3, double& p4, NumericVector& expt, NumericVector& s_i, NumericVector& info, int se, int bnb);
RcppExport SEXP _bzinb_dBvZINB_Expt_vec(SEXP xvecSEXP, SEXP yvecSEXP, SEXP freqSEXP, SEXP nSEXP, SEXP a0SEXP, SEXP a1SEXP, SEXP a2SEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP p1SEXP, SEXP p2SEXP, SEXP p3SEXP, SEXP p4SEXP, SEXP exptSEXP, SEXP s_iSEXP, SEXP infoSEXP, SEXP seSEXP, SEXP bnbSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type xvec(xvecSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type yvec(yvecSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type freq(freqSEXP);
    Rcpp::traits::input_parameter< int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< double& >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< double& >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< double& >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< double& >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< double& >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< double& >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< double& >::type p2(p2SEXP);
    Rcpp::traits::input_parameter< double& >::type p3(p3SEXP);
    Rcpp::traits::input_parameter< double& >::type p4(p4SEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type expt(exptSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type s_i(s_iSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type info(infoSEXP);
    Rcpp::traits::input_parameter< int >::type se(seSEXP);
    Rcpp::traits::input_parameter< int >::type bnb(bnbSEXP);
    dBvZINB_Expt_vec(xvec, yvec, freq, n, a0, a1, a2, b1, b2, p1, p2, p3, p4, expt, s_i, info, se, bnb);
    return R_NilValue;
END_RCPP
}
// inv_digamma
double inv_digamma(double x, double y);
RcppExport SEXP _bzinb_inv_digamma(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(inv_digamma(x, y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bzinb_em", (DL_FUNC) &_bzinb_em, 10},
    {"_bzinb_dBvZINB_Expt", (DL_FUNC) &_bzinb_dBvZINB_Expt, 17},
    {"_bzinb_dBvZINB_Expt_vec", (DL_FUNC) &_bzinb_dBvZINB_Expt_vec, 18},
    {"_bzinb_inv_digamma", (DL_FUNC) &_bzinb_inv_digamma, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_bzinb(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
