/**
 * $Id$
 * Author: Stefano Carrazza, stefano.carrazza@mi.infn.it
 */

#pragma once

#include <vector>
#include <string>

using namespace std;

enum {
  ER_NONE,   //!< PDF set is not an error set
  ER_MC,     //!< 1sigma error for NNPDF
  ER_EIG,    //!< 1sigma error for CTEQ & MSTW
  ER_EIG90,  //!< 90cl error for CTEQ & MSTW
  ER_SYMEIG  //!< 68cl sym.eigenstate ABM
};

double ComputeAVG(int n, const double *x);             //!< Compute average from x points
double ComputeAVG(int n, int pdf, const double *x);    //!< Compute average from x points
double ComputeAVG(vector<double> x);                 //!< Compute average from vector<double> x
double ComputeStdDev(int n, const double *x);          //!< Compute the std deviation
double ComputeStdDev(int n, int pdf, const double *x); //!< Compute the std deviation
double ComputeStdDev(vector<double> x); //!< Compute the std deviation
double ComputeEigErr(int p, const double *x);          //!< Compute error in the Hessian method for PDFs
double ComputeEigErr(int p, int n,const double *x);    //!< Compute error in the Hessian method for Obs
double ComputeSymEigErr(int p, double cv, const double *x);          //!< Compute error in the Hessian method for PDFs
