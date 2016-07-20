/**
 * Fortran WRAPPER FOR HELL
 * stefano.carrazza@mi.infn.it
 */

#include "HELL/include/hell.hh"
#include "HELL/include/hell_wrapper.hh"

#include <iostream>
#include <string>
using namespace std;

#define STR_EXPAND(top) #top
#define STR(tok) STR_EXPAND(tok)

// Global allocation
HELLx::LogOrder order;
HELLx::HELLx *sxD = NULL;
HELLx::Order fixed_order_to_be_matched_to;
HELLx::sqmatrix<double> xdPNLL;

string HELLdataPath()
{
  string dataDir(STR(DATA_PATH));
  stringstream s;
  s << dataDir << "/apfel";
  return s.str();
}

extern "C" {

  void helllogorder_(int *ord)
  {    
    if (*ord == 0)
      order = HELLx::LL;
    else 
      order = HELLx::NLL;
  }

  void hell_()
  {
    if (sxD) delete sxD;
    sxD = new HELLx::HELLx(order, HELLdataPath());
  }

  void hellorder_(int *ord)
  {
    if (*ord == -1)
      fixed_order_to_be_matched_to = HELLx::none;
    else if (*ord == 0)
      fixed_order_to_be_matched_to = HELLx::LO;
    else if (*ord == 1)
      fixed_order_to_be_matched_to = HELLx::NLO;
    else if (*ord == 2)
      fixed_order_to_be_matched_to = HELLx::NNLO;	
  }

  // Splitting functions
  double xdeltap_(int *nf, int *k, double *as, double *x)
  {
    xdPNLL = sxD->DeltaP(*nf, *as, *x, fixed_order_to_be_matched_to);
    double xdp = 0;
    if(*k == 4)      xdp = xdPNLL.qq();
    else if(*k == 5) xdp = xdPNLL.qg();
    else if(*k == 6) xdp = xdPNLL.gq();
    else if(*k == 7) xdp = xdPNLL.gg();
    return xdp;
  }

  // Coeffincient functions
  double xdeltac2_(int *nf, int *k, double *as, double *x)
  {
    double xdc2 = 0;
    if(*k == 1)      xdc2 = sxD->deltaC2g(*nf, *as, *x, fixed_order_to_be_matched_to); // Gluon
    else if(*k == 2) xdc2 = sxD->deltaC2q(*nf, *as, *x, fixed_order_to_be_matched_to); // Quark
    return xdc2;
  }

  double xdeltacl_(int *nf, int *k, double *as, double *x)
  {
    double xdcl = 0;
    if(*k == 1)      xdcl = sxD->deltaCLg(*nf, *as, *x, fixed_order_to_be_matched_to); // Gluon
    else if(*k == 2) xdcl = sxD->deltaCLq(*nf, *as, *x, fixed_order_to_be_matched_to); // Quark
    return xdcl;
  }

}
