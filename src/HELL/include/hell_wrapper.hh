/**
 * Fortran WRAPPER FOR HELL
 * stefano.carrazza@mi.infn.it
 */

#pragma once

#ifndef DATA_PATH
#define DATA_PATH ./data
#endif

extern "C" {
  // Initializers
  void helllogorder_(int *ord);
  void hell_();
  void hellorder_(int *ord);

  // Splitting functions
  double xdeltap_(int *nf, int *k, double *as, double *x);

  // Matching coefficients
  double xdeltak_ (int *nf, int *k, double *as, double *x, double *m_Q_ratio);

  // Massless coefficient functions
  double xdeltac2_(int *nf, int *k, double *as, double *x);
  double xdeltacl_(int *nf, int *k, double *as, double *x);

  // Massive coefficient functions
  double xdeltamc2_(int *nf, int *k, double *as, double *x, double *m_Q_ratio);
  double xdeltamcl_(int *nf, int *k, double *as, double *x, double *m_Q_ratio);

}

