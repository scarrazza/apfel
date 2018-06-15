#ifndef APFELEVOL_H
#define APFELEVOL_H

#include <string>
#include <iostream>
using std::string;

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

  /// Computes evolution 
  void EvolveAPFEL(double Q0, double Q);

  /// Computes the derivative of PDFs 
  void DeriveAPFEL(double Q);

  /// Cache PDFs on a (x,Q2)-grid
  void CachePDFsAPFEL(double Q0);

  /// Returns x*PDF
  double xPDF(int i, double x);

  /// Returns x*PDF has functions of x and Q using the cached evolution
  double xPDFxQ(int i, double x, double Q);

  /// Returns the derivative of x*PDF
  double dxPDF(int i, double x);

  /// Returns x*PDF on the joint grid
  double xPDFj(int i, double x);

  /// Returns x*gamma
  double xgamma(double x);

  /// Returns x*gamma on the joint grid
  double xgammaj(double x);

  /// Returns the derivative of x*gamma
  double dxgamma(double x);

  /// Returns all PDFs on the joint grid
  void xPDFall(double x, double *xf);

  /// Returns all PDFs (including the photon) on the joint grid
  void xPDFallPhoton(double x, double *xf);

  /// Returns all the cached PDFs on the (x,Q)-grid
  void xPDFxQall(double x, double Q, double *xf);

  /// Returns x*Lepton PDFs
  double xLepton(int i, double x);

  /// Returns x*Lepton PDFs on the joint grid
  double xLeptonj(int i, double x);

  /// External Evolution Operator
  double ExternalEvolutionOperator(const std::string& fname, int i, int j, double x, int beta);

  /// External Evolution Matrix evolution to evolution
  double ExternalEvolutionMatrixEv2Ev(int i, int j, int alpha, int beta);

  /// External Evolution Matrix evolution to physical
  double ExternalEvolutionMatrixEv2Ph(int i, int j, int alpha, int beta);

  /// External Evolution Matrix evolution to physical
  double ExternalEvolutionMatrixPh2Ph(int i, int j, int alpha, int beta);

  /// Compute External Splitting Functions
  void ComputeExternalSplittingFunctions(const std::string& fname, int pt, int nf, double x, int beta);

  /// Return External Splitting Functions
  double ExternalSplittingFunctions(int i, int j);

  /// Builds the *.LHgrid output file
  void LHAPDFgrid(int Nrep, double Qin, const std::string& fname);

  /// Builds the *.LHgrid output file with the derivative of the input set
  void LHAPDFgridDerivative(int Nrep, const std::string& fname);

  /// Returns the value of alpha_QCD at the given scale
  double AlphaQCD(double Q);

  /// Returns the value of alpha_QED at the given scale
  double AlphaQED(double Q);

  /// Returns Heavy Quark Masses
  double HeavyQuarkMass(int,double);

  /// Returns Heavy Quark Threholds
  double GetThreshold(int);

  /// Gets the maximum number of flavours that the evolution 
  /// of alphaQCD and alphaQED can reach
  int GetMaxFlavourAlpha();

  /// Gets the maximum number of flavours that the evolution 
  /// of PDFs can reach.
  int GetMaxFlavourPDFs();

  /// Returns Heavy Quark Threholds (identical to GetThreshold)
  double HeavyQuarkThreshold(int);

  /// Returns the N-th Mellin moment of the i-th PDF 
  /// in the physical basis at the final scale
  double NPDF(int i, int N);

  /// Returns the N-th Mellin moment of the photon PDF 
  /// in the physical basis at the final scale
  double Ngamma(int N);

  /// Returns the Luminosity 
  double LUMI(int i, int j, double S);

  /// Returns the joint x-space grid 
  double xGrid(int alpha);

  /// Returns the number of intervals of the joint x-space grid
  int nIntervals();

  /// Gets APFEL version
  std::string GetVersion(void);

  /*
   * Set Methods
   */

  /// Cleans up parameters
  void CleanUp(void);

  /// Enables welcome message
  void EnableWelcomeMessage(int);

  /// Enables evolution operator computation
  void EnableEvolutionOperator(int);

  /// Enables the evolution of leptons
  void EnableLeptonEvolution(int);

  /// Lock internal subgrids
  void LockGrids(int);

  /// Switch to the time-like evolution
  void SetTimeLikeEvolution(int);

  /// Switch to the polarized evolution
  void SetPolarizedEvolution(int);

  /// Switch to the fast evolution
  void SetFastEvolution(int);

  /// Enables the running of the MSbar masses
  void EnableMassRunning(int);

  /// Switch on the small-x resummation
  void SetSmallxResummation(int, const std::string& la);  

  /// Sets the reference values of $alpha_{s}$ at the reference scale
  void SetAlphaQCDRef(double alpharef, double Qref);

  /// Sets the reference values of $alpha$ at the reference scale
  void SetAlphaQEDRef(double alpharef, double Qref);

  /// Sets the solution of the beta function
  void SetAlphaEvolution(const std::string& evol);

  /// Sets the value of LambdaQCD for "nref" flavours
  void SetLambdaQCDRef(double lambdaref, int nref);

  /// Sets the solution of the DGLAP equation
  void SetPDFEvolution(const std::string& evolp);

  /// Sets the epsilon parameter of the truncated DGLAP solution
  void SetEpsilonTruncation(double eps);

  /// Sets the minimimum and the maximum energy allowed for the evolution
  void SetQLimits(double Qmin, double Qmax);

  /// Sets the FFNS as a default
  void SetFFNS(int nfl);

  /// Sets the parameters of the i-th x-space grid
  void SetGridParameters(int i, int np, int deg, double x);

  /// Sets the parameters of the Q-space grid
  void SetQGridParameters(int npq, int degq);

  /// Sets the parameters of the grid over which PDFs in the LHAPDF
  /// format will be witten
  void SetLHgridParameters(int nx, int nxm, double xmin, double xm, double xmax,
			   int nq2, double q2min, double q2max);

  /// Sets the user given i-th x-space grid
  void SetExternalGrid(int i, int np, int deg, double *x);

  /// Sets the maximum number of flavours that the evolution 
  /// of alphaQCD and alphaQED can reach
  void SetMaxFlavourAlpha(int nf);

  /// Sets the maximum number of flavours that the evolution 
  /// of PDFs can reach.
  void SetMaxFlavourPDFs(int nf);

  /// Sets as a default the heavy quark MSbar masses
  void SetMSbarMasses(double mc, double mb, double mt);

  /// Sets scales at which heavy quark masses are defined.
  /// No effect for the pole masses.
  void SetMassScaleReference(double Qc, double Qb, double Qt);

  /// Sets ratios between heavy quark masses and thresholds 
  void SetMassMatchingScales(double kmc, double kmb, double kmt);
  
  /// Sets the number of x-space grids that will be used in the computation
  void SetNumberOfGrids(int n);

  /// Sets the name of the PDF set to be used at the initial scale
  void SetPDFSet(const std::string& name);

  /// Sets the perturbative order of the evolution
  void SetPerturbativeOrder(int pto);

  /// Returns the perturbative order of the evolution
  int GetPerturbativeOrder();

  /// Returns the final factorization scale
  double GetMuF();

  /// Returns the initial factorization scale
  double GetMuF0();

  /// Sets as a default the heavy quark pole masses
  void SetPoleMasses(double mc, double mb, double mt);

  /// Set the value of the tau mass in GeV
  void SetTauMass(double masst);

  /// Sets the ratio between renormalization and factorization scales.
  void SetRenFacRatio(double ratio);
  
  /// Sets the replica to be used as initial PDFs (only with a LHAPDF grid)
  void SetReplica(int nr);
  
  /// Sets the FFNS as a default
  void SetTheory(const std::string& theory);
  
  /// Enables NLO QED corrections
  void EnableNLOQEDCorrections(int);
  
  /// Sets the VFNS as a default
  void SetVFNS(void);

  /// List the functions of APFEL
  void ListFunctions(void);

  /// Check APFEL
  bool CheckAPFEL(void);
}

#endif
