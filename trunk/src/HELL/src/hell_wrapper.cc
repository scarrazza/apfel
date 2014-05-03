/**
 * Fortran WRAPPER FOR HELL
 * stefano.carrazza@mi.infn.it
 */

#include "HELL/include/hell.hh"
#include "HELL/include/hell_wrapper.hh"

#include <iostream>
#include <string>
using namespace std;


// Global allocation
HELL::LogOrder order;
HELL::HELL *sxD = NULL;
HELL::Order fixed_order_to_be_matched_to;
sqmatrix<double> xdPNLL;

extern "C" {

  void helllogorder_(int *ord)
  {    
    if (*ord == 0)
      order = HELL::LL;
    else 
      order = HELL::NLL;
  }

  void hell_(double *asmc, double *asmb, double *asmt)
  {
    if (sxD) delete sxD;
    sxD = new HELL::HELL(order, "./data/");
    sxD->init_as_thresholds(*asmc,*asmb,*asmt);
  }

  void hellorder_(int *ord)
  {
    if (*ord == -1)
      fixed_order_to_be_matched_to = HELL::none;
    else if (*ord == 0)
      fixed_order_to_be_matched_to = HELL::LO;
    else if (*ord == 1)
      fixed_order_to_be_matched_to = HELL::NLO;
    else if (*ord == 2)
      fixed_order_to_be_matched_to = HELL::NNLO;	
  }

  void xdeltap_(double *as, double *x)
  {
    xdPNLL = sxD->xDeltaP(*as, *x, fixed_order_to_be_matched_to);
    cout << endl << "Printing x*DeltaP_ij(x=" << *x << ", as=" << *as << ")"  << endl;
    cout << "gg: " << xdPNLL.gg() << endl
	 << "gq: " << xdPNLL.gq() << endl
	 << "qg: " << xdPNLL.qg() << endl
	 << "qq: " << xdPNLL.qq() << endl;
  }
}
