/**
 * Fortran WRAPPER FOR HELL
 * stefano.carrazza@mi.infn.it
 */

#pragma once

#ifndef DATA_PATH
#define DATA_PATH ./data
#endif

extern "C" {
  void helllogorder_(int *ord);
  void hell_(double *asmc, double *asmb, double *asmt);
  void hellorder_(int *ord);
  void xdeltap_(double *as, double *x);
}

