#ifndef APFELOBS_H
#define APFELOBS_H

#include <string>
#include <iostream>

#include "APFEL/APFELevol.h"

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
   * DIS module
   */
  /// Initialize the new DIS module
  void InitializeAPFEL_DIS();

  /// Precompute the structure functions
  void ComputeStructureFunctionsAPFEL(double Q0, double Q);

  /// Cache structure functions on a (x,Q2)-grid
  void CacheStructureFunctionsAPFEL(double Q0);

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

  /// Select a given charge contribution
  void SelectCharge(const std::string& selch);

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

  /// Polarised structure functions
  double g1light(double x);
  double g1charm(double x);
  double g1bottom(double x);
  double g1top(double x);
  double g1total(double x);

  double gLlight(double x);
  double gLcharm(double x);
  double gLbottom(double x);
  double gLtop(double x);
  double gLtotal(double x);

  double g4light(double x);
  double g4charm(double x);
  double g4bottom(double x);
  double g4top(double x);
  double g4total(double x);

  /// Chached structure functions
  double StructureFunctionxQ(const std::string& proc,const std::string& sf,const std::string& comp,double x,double Q);

  /// Set the value of the Z mass in GeV
  void SetZMass(double massz);

  /// Set the value of the W mass in GeV
  void SetWMass(double massw);

  /// Set the value of the proton mass in GeV
  void SetProtonMass(double massp);

  /// Set the value of sin^2(theta_W)
  void SetSin2ThetaW(double sw);

  /// Set the absolute value of the entries of the CKM matrix
  void SetCKM(double vud,double vus,double vub,double vcd,double vcs,double vcb,double vtd,double vts,double vtb);

  /// Set the the correction to the Z propagator
  void SetPropagatorCorrection(double dr);

  /// Set the EW vector and axial couplings
  void SetEWCouplings(double vd,double vu,double ad,double au);

  /// Set the value of the Fermi constant
  void SetGFermi(double gf);

  /// Sets the ratio between renormalization scale and Q
  void SetRenQRatio(double ratioR);

  /// Sets the ratio between factorization scale and Q
  void SetFacQRatio(double ratioF);

  /// Enables the possibility to vary scales dynamically
  void EnableDynamicalScaleVariations(int);

  /// Enables intrinsic charm contributions
  void EnableIntrinsicCharm(int);

  // Returns the mass of the Z
  double GetZMass();

  // Returns the mass of the W
  double GetWMass();

  // Returns the mass of the proton
  double GetProtonMass();

  // Returns sin^2(\theta_W)
  double GetSin2ThetaW();

  // Returns the entries of the CKM matrix
  double GetCKM(int u, int d);

  // Returns G_{Fermi}
  double GetGFermi();

  // Returns the SIA total cross section
  double GetSIATotalCrossSection(int pto, double q, const std::string& comp);

  /// Enables the target mass corrections
  void EnableTargetMassCorrections(int);

  /// Enables the FONLL damping factor
  void EnableDampingFONLL(int);

  /// Set the FONLL damping power suppression
  void SetDampingPowerFONLL(int, int, int);

  /// Get the EW DIS charges at the scale Q2
  void ComputeChargesDIS(double q2, double *bq, double *dq, double *bqt);

  /// Get the EW DIS charges at the scale Q2
  void ComputeChargesDIS(double q2, double *bq, double *dq, double *bqt);

  /// Compute F2 at LO using PDFs
  double F2LO(double x, double q);

  /// Emulator of the FKgenerator
  double FKSimulator(double x,double q,double y,int i,int beta);

  /// Set the observable for FKgenerator simulator
  void SetFKObservable(const std::string& obs);

  /// Get the observable used by the FKgenerator simulator
  void GetFKObservable();

  /// Observable according to the FKgenerator naming
  double FKObservables(double x,double q,double y);
						
  /// Functions for FTDY
  void ComputeFKTables(const std::string& inputfile, const std::string& outputpath,
		       double Q0, int* flmap);

  void ComputeHardCrossSectionsDY(const std::string& datafile, 
				  const std::string& outputfile);

  /// Enables NLO QED corrections in the structure functions
  void EnableSFNLOQEDCorrections(int);

  /// Builds the *.LHgrid output file
  void LHAPDFgridStructureFunctions(int Nrep, double Qin, const std::string& fname);

  /// Set procedure for the scale variation
  void SetScaleVariationProcedure(int svp);
}

#endif
