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
  void hell_();
  void hellorder_(int *ord);
  double xdeltap_(int *nf, int *k, double *as, double *x);
  double xdeltac2_(int *nf, int *k, double *as, double *x);
  double xdeltacl_(int *nf, int *k, double *as, double *x);
}

