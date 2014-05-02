/**
 * Fortran WRAPPER FOR HELL
 * stefano.carrazza@mi.infn.it
 */

#pragma once
#include "HELL/include/hell.hh"

extern "C" {

  void HELLLogOrder_(int *ord);
  void HELL_(double *asmc, double *asmb, double *asmt);
  void HELLOrder_(int *ord);
  void xDeltaP_(double *as, double *x);

}
