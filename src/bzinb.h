#ifndef BZINB_H
#define BZINB_H

#include <Rcpp.h>
#include <math.h>
#include <string>
#include <iostream>
#include <numeric>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
// 
// double EPSILON1;
// double EPSILON2;

void l1(int& x, int& y, long double& a0, long double& a1, long double& a2, int &k, int& m, 
        long double& result, long double adjj = 0);
void l1_c (long double& t1, long double& t2, int& k, int& m, long double& result, long double adjj);
void l1_AC (long double& t1, long double& t2, int& x, int& y, long double& a0, long double& a1, 
            long double& a2, int& k, int& m, long double& result, long double adjj = 0);
void l2_A (int& x, long double& a0, long double& a1, long double& a2, int& k, 
           long double& result, long double adjj);
void l3_A (int& y, long double& a0, long double& a1, long double& a2, int& m, 
           long double& result, long double adjj);
void R0_E1(int& x, int& y, int& k, int& m, long double& a0, long double& result);
long double log_R0_E1(int& x, int& y, int& k, int& m, long double& a0);
long double	log_R0_E2(int& x, long double& a0, int& k);
long double	log_R0_E3(int& y, long double& a0, int& m);
void R1_E1(int& k, long double& a1, long double& result);
long double log_R1_E1(int& k, long double& a1);

void dBvZINB_Expt(int &x, int &y, int &freq, long double &a0, long double &a1, long double &a2,
                  long double &b1, long double &b2, long double &p1, long double &p2, 
                  long double &p3, long double &p4,
                  Rcpp::NumericVector &expt, Rcpp::NumericVector &s_i, Rcpp::NumericVector &info, int se);
void dBvZINB_Expt_vec(const Rcpp::IntegerVector &xvec, const Rcpp::IntegerVector &yvec, 
                      const Rcpp::IntegerVector &freq, 
                      const int &n, long double &a0, long double &a1, long double &a2,
                      long double &b1, long double &b2, long double &p1, long double &p2, 
                      long double &p3, long double &p4,
                      Rcpp::NumericVector &expt, Rcpp::NumericVector &s_i, Rcpp::NumericVector &info, int se);

long double inv_digamma(long double x, long double y);
void inv_digamma_vec(long double lb[1], Rcpp::NumericVector &expt, long double a[3], 
                     long double idgam[3]);
long double hfunc(long double lb[1], Rcpp::NumericVector &expt, long double a[3], 
                  long double idgam[3]);
void opt_lb(long double lb[1], Rcpp::NumericVector &expt, long double a[3], 
            long double idgam[3]);

void em(Rcpp::NumericVector& param2, const Rcpp::IntegerVector &xvec, const Rcpp::IntegerVector &yvec, 
        const Rcpp::IntegerVector &freq, const int &n, Rcpp::NumericVector &expt, Rcpp::NumericVector &info,
        const int &se, Rcpp::IntegerVector &iter, int &maxiter, double &tol, int showFlag, 
        Rcpp::IntegerVector &nonconv, Rcpp::NumericVector trajectory);

#endif
