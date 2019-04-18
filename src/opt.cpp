// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>
#include "bzinb.h"
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
using namespace std;
using namespace Rcpp;

#define EPSILON1 1e-9  //for inverse-gamma
#define EPSILON2 1e-9  //for M-step
// #define DEBUG
// #define DEBUG2
// #define DEBUG3
//bool debug = true;

// 2. optimization
// [[Rcpp::export]]
double inv_digamma(double x, double y) 
{ 
  double h = (boost::math::digamma(x) - y) / boost::math::trigamma(x);
  while (fabs(h) >= EPSILON1)
  {
    h = (boost::math::digamma(x) - y) / boost::math::trigamma(x);
    while (h > x) {
      #ifdef DEBUG 
        Rcout <<"inv_digamm hit zero:" << h << " ";
      #endif
      h /= 2.0;
    }
    x = x - h;
  }
  // if (debug) {Rcout << endl << "inv_digamm = " << x << endl;}
  #ifdef DEBUG 
    Rcout << "inv_digamm = " << x << " ";
  #endif
  
  return(x);
}


// optimization


void inv_digamma_vec(double lb[1], NumericVector &expt, double a[3], double idgam[3]) 
{ 
  //double idgam[3];
  //double *idgam = new double[3];
  idgam[0] = inv_digamma(a[0], expt[4] - lb[0]);
  idgam[1] = inv_digamma(a[1], expt[5] - lb[0]);
  idgam[2] = inv_digamma(a[2], expt[6] - lb[0]);
  //return(idgam);
}

// double objective(double lb, NumericVector &expt, NumericVector &a, NumericVector &idgam) 
// { 
//   inv_digamma_vec(lb, expt, a, idgam);
//   // Rcout << "idgam123: "<< idgam[0] << " "<< idgam[1] << " "<< idgam[2] << endl;
//   //  inv_digamma_vec(b1, expt, a, idgam);
//   // double result = (double);
//   return (log(idgam[0] + idgam[1] + idgam[2]) + lb - log(expt[1] + expt[2] + expt[3]));
// }
// 
// // Derivative
// double derivFunc(double lb, NumericVector &expt, NumericVector &a, NumericVector &idgam)
// {
//   inv_digamma_vec(lb, expt, a, idgam);
//   double result = 0.0;
//   for (int i = 0; i < 3; i++) {
//     result += 1/ boost::math::trigamma(idgam[i]);
//   }
//   return ((result)/(idgam[0] + idgam[1] + idgam[2]) + 1.0);
// }

// Derivative
double hfunc(double lb[1], NumericVector &expt, double a[3], double idgam[3])
{
  inv_digamma_vec(lb, expt, a, idgam);
  double result = 0.0;
  double obj = (log(idgam[0] + idgam[1] + idgam[2]) + lb[0] - log(expt[1] + expt[2] + expt[3]));
#ifdef DEBUG 
  Rcout << "objective = " << obj << " ";
#endif
  for (int i = 0; i < 3; i++) {
    // Rcout << "idgam[0:2]" << idgam[0] << " " << idgam[1] << " " << idgam[2] << endl;
    result += (1/ boost::math::trigamma(idgam[i]));
  }
  result = - result / (idgam[0] + idgam[1] + idgam[2]) + 1.0;
  result = obj / result;
  return (result);
}

void opt_lb(double lb[1], NumericVector &expt, double a[3], double idgam[3])
{
  int i = 1;
  //double* lb = log(b1);
  // double h =  objective(lb, expt, a, idgam) / derivFunc(lb, expt, a, idgam);
  double h = hfunc(lb, expt, a, idgam);
  // if (debug) {Rcout << "h = " << h << " ";}
#ifdef DEBUG 
  Rcout << "h = " << h << " ";
#endif
#ifdef DEBUG3
  Rcout << "opt_iter = 1 ";  
#endif
  while ((fabs(h) >= EPSILON2) & (i < 11))  //sometimes h oscilliates around zero 
                // (usually close to solution with only a few iterations, as the initial value is close to truth.)
                // See p.10 of http://mathfaculty.fullerton.edu/mathews//n2003/newtonsmethod/Newton'sMethodProof.pdf
  {
    i += 1;
    // h = objective(lb, expt, a, idgam)/derivFunc(lb, expt, a, idgam);
    // while (h > b1) {
    //   // Rcout << "opt_b1 hit zero. ";
    //   h /= 2.0;
    // }
    // Rcout << "h = " << h << " ";
#ifdef DEBUG3 
    Rcout << "(i= " << i<< ", lb= " << lb[0] << ", h= " << h << ", lb-h=" << lb[0]-h << ").. ";
#endif
    lb[0] -= h;
    h = hfunc(lb, expt, a, idgam);
  }
#ifdef DEBUG2 
  Rcout << "lb =" << lb[0]  << ", b1 = " << exp(lb[0]) << " / ";
#endif
  
// if (debug) {Rcout << "lb =" << lb  << ", b1 = " << exp(lb) << endl;}
// if (debug) {Rcout << "lb =" << lb  << endl;}
  // double* idgam = inv_digamma_vec(lb, expt, a);
  a[0] = idgam[0];
  a[1] = idgam[1];
  a[2] = idgam[2];
  // Rcout << "a123: "<< a[0] << " "<< a[1] << " "<< a[2] << " " << exp(lb) << endl;
}
