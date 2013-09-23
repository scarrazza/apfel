#ifndef APFEL_H
#define APFEL_H

#include <string>
#include <iostream>

/**
 * @mainpage A C++ wrapper for the APFEL library
 *
 * @section intro Introduction
 * The APFEL library provides a set of C++ wrapper functions
 * for its Fortran subroutines.
 */

/// Namespace containing all the APFEL wrapper functions.
namespace APFEL {

  /*
   * Main Methods
   */ 

  /// Initialize the library
  void InitializeAPFEL(void);

  /// Compute evolution 
  void EvolveAPFEL(double Q0, double Q);

  /// Return x*PDF
  double xPDF(int i, double x);

  /// Return x*gamma
  double xgamma(double x);

  /// Build the *.LHgrid output file
  void LHAPDFgrid(int Nrep, double Qin, const std::string& fname);

  /// Returns the value of alpha_QCD at the given scale
  double AlphaQCD(double Q);

  /// Returns the value of alpha_QED at the given scale
  double AlphaQED(double Q);

  /// Returns the N-th Mellin moment of the i-th PDF 
  /// in the physical basis at the final scale
  double NPDF(int i, int N);

  /// Returns the N-th Mellin moment of the photon PDF 
  /// in the physical basis at the final scale
  double Ngamma(int N);
  
  /*
   * Set Methods
   */

  /// Set the reference values of $alpha_{s}$ at the reference scale
  void SetAlphaQCDRef(double alpharef, double Qref);

  /// Set the reference values of $alpha$ at the reference scale
  void SetAlphaQEDRef(double alpharef, double Qref);

  /// Set the minimimum and the maximum energy allowed for the evolution
  void SetQLimits(double Qmin, double Qmax);

  /// Set the FFNS as a default
  void SetFFNS(int nfl);

  /// Set the parameters of the i-th x-space grid
  void SetGridParameters(int i, int np, int deg, double x);

  /// Set the maximum number of flavours that the evolution 
  /// of alphaQCD and alphaQED can reach
  void SetMaxFlavourAlpha(int nf);

  /// Set the maximum number of flavours that the evolution 
  /// of PDFs can reach.
  void SetMaxFlavourPDFs(int nf);

  /// Set as a default the heavy quark MSbar masses
  void SetMSbarMasses(double mc, double mb, double mt);
  
  /// Set the number of x-space grids that will be used in the computation
  void SetNumberOfGrids(int n);

  /// Set the name of the PDF set to be used at the initial scale
  void SetPDFSet(const std::string& name);

  /// Set the perturbative order of the evolution
  void SetPerturbativeOrder(int pto);

  /// Set as a default the heavy quark pole masses
  void SetPoleMasses(double mc, double mb, double mt);

  /// Set the ratio between renormalization and factorization scales.
  void SetRenFacRatio(double ratio);
  
  /// Set the replica to be used as initial PDFs (only with a LHAPDF grid)
  void SetReplica(int nr);
  
  /// Set the FFNS as a default
  void SetTheory(const std::string& theory);
  
  /// Set the VFNS as a default
  void SetVFNS(void);

}

#endif
