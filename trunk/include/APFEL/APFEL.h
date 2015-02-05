#ifndef APFEL_H
#define APFEL_H

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

  /// Returns x*PDF
  double xPDF(int i, double x);

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

  /// External Evolution Operator
  double ExternalEvolutionOperator(const std::string& fname, int i, int j, double x, int beta);

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

  /// Lock internal subgrids
  void LockGrids(int);

  /// Switch to the time-like evolution
  void SetTimeLikeEvolution(int);

  /// Switch to the fast evolution
  void SetFastEvolution(int);

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

  /// Sets the minimimum and the maximum energy allowed for the evolution
  void SetQLimits(double Qmin, double Qmax);

  /// Sets the FFNS as a default
  void SetFFNS(int nfl);

  /// Sets the parameters of the i-th x-space grid
  void SetGridParameters(int i, int np, int deg, double x);

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
  
  /// Sets the number of x-space grids that will be used in the computation
  void SetNumberOfGrids(int n);

  /// Sets the name of the PDF set to be used at the initial scale
  void SetPDFSet(const std::string& name);

  /// Sets the perturbative order of the evolution
  void SetPerturbativeOrder(int pto);

  /// Sets as a default the heavy quark pole masses
  void SetPoleMasses(double mc, double mb, double mt);

  /// Sets the ratio between renormalization and factorization scales.
  void SetRenFacRatio(double ratio);
  
  /// Sets the replica to be used as initial PDFs (only with a LHAPDF grid)
  void SetReplica(int nr);
  
  /// Sets the FFNS as a default
  void SetTheory(const std::string& theory);
  
  /// Sets the VFNS as a default
  void SetVFNS(void);

  /*
   * DIS module
   */

  /// DIS observables
  void DIS_xsec(double x,double qi,double qf,double y,double pol,
		const std::string& proc,const std::string& scheme,
		int pto,const std::string& pdfset, int irep,
		const std::string& target, const std::string& proj,
		double *F2, double *F3, double *FL, double *sigma);

  /*
   * New DIS module
   */

  /// Initialize the new DIS module
  void InitializeAPFEL_DIS();

  /// Precompute the structure functions
  void ComputeStructureFunctionsAPFEL(double Q0, double Q);

  /// Set the mass scheme for the structure functions
  void SetMassScheme(const std::string& ms);

  /// Set the polarization
  void SetPolarizationDIS(double pol);

  // Set the process of the structure functions (EM, NC or CC)
  void SetProcessDIS(const std::string& pr);

  /// Set the projectile
  void SetProjectileDIS(const std::string& lept);

  /// Set the target
  void SetTargetDIS(const std::string& tar);

  /// Returns the DIS operator times the evolution factors on the grid
  double ExternalDISOperator(const std::string& SF,int ihq,int i,double x,int beta);

  /// Structure functions
  double F2light(double x);
  double F2charm(double x);
  double F2bottom(double x);
  double F2top(double x);
  double F2total(double x);

  double FLlight(double x);
  double FLcharm(double x);
  double FLbottom(double x);
  double FLtop(double x);
  double FLtotal(double x);

  double F3light(double x);
  double F3charm(double x);
  double F3bottom(double x);
  double F3top(double x);
  double F3total(double x);

  /// Set the value of the Z mass in GeV
  void SetZMass(double massz);

  /// Set the value of the W mass in GeV
  void SetWMass(double massw);

  /// Set the value of the proton mass in GeV
  void SetProtonMass(double massp);

  /// Set the value of sin(theta_W)
  void SetDinThetaW(double sw);

  /// Set the absolute value of the entries of the CKM matrix
  void SetCKM(double vud,double vus,double vub,double vcd,double vcs,double vcb,double vtd,double vts,double vtb);

  /// Set the value of the Fermi constant
  void SetGFermi(double gf);

  /// Emulator of the FKgenerator
  double FKSimulator(const std::string& obs,double x,double q,double y,int i,int beta);

}

#endif
