/*
  TabulationFragCxx.cc:

  Example program used for fragmentation function evolution.
  The fragmentation functions are hardcoded and for the following
  options are available:
                          
  DSS                  HKNS:
                          
  - "dss_pip_lo"       - "hknsff07_pip_lo"
  - "dss_pip_nlo"      - "hknsff07_pip_nlo"
  - "dss_pim_lo"       - "hknsff07_pim_lo"
  - "dss_pim_nlo"      - "hknsff07_pim_nlo"
  - "dss_Kp_lo"        - "hknsff07_Kp_lo"
  - "dss_Kp_nlo"       - "hknsff07_Kp_nlo"
  - "dss_Km_lo"        - "hknsff07_Km_lo"
  - "dss_Km_nlo"       - "hknsff07_Km_nlo"
  - "dss_p_lo"         - "hknsff07_p_lo"
  - "dss_p_nlo"        - "hknsff07_p_nlo"
  - "dss_pb_lo"        - "hknsff07_pb_lo"
  - "dss_pb_nlo"       - "hknsff07_pb_nlo"
  - "dss_h_lo"    
  - "dss_h_nlo"   
  - "dss_hb_lo"   
  - "dss_hb_nlo"  
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "APFEL/APFEL.h"
using namespace std;

int main()
{
  double xlha[] = {1e-2, 5e-2, 1e-1, 2e-1, 3e-1, 4e-1, 
		5e-1, 6e-1, 7e-1, 8e-1, 9e-1};

  //
  // Settings
  //
  APFEL::SetTimeLikeEvolution(true);
  //
  // HKNS
  //
  APFEL::SetPDFSet("hknsff07_p_lo");
  APFEL::SetPerturbativeOrder(0);
  APFEL::SetMaxFlavourPDFs(5);
  APFEL::SetMaxFlavourAlpha(5);
  APFEL::SetPoleMasses(1.43,4.3,175);
  APFEL::SetAlphaEvolution("lambda");
  APFEL::SetLambdaQCDRef(0.220,4);
  //APFEL::SetLambdaQCDRef(0.323,4);
  APFEL::SetNumberOfGrids(2);
  APFEL::SetGridParameters(1,50,3,1e-2);
  APFEL::SetGridParameters(2,40,3,7e-1);
  /*
  //
  // DSS
  //
  APFEL::SetPDFSet("dss_p_lo");
  APFEL::SetPerturbativeOrder(0);
  APFEL::SetMaxFlavourPDFs(5);
  APFEL::SetMaxFlavourAlpha(5);
  APFEL::SetPoleMasses(1.43,4.3,175);
  APFEL::SetAlphaEvolution("lambda");
  APFEL::SetLambdaQCDRef(0.220,4);
  //APFEL::SetLambdaQCDRef(0.334,4);
  APFEL::SetNumberOfGrids(2);
  APFEL::SetGridParameters(1,50,3,5e-2);
  APFEL::SetGridParameters(2,40,3,7e-1);
  */
  // Initializes integrals on the grids
  APFEL::InitializeAPFEL();

  double Q02, Q2, eps = 1e-10;
  cout << "Enter initial and final scale in GeV^2" << endl;
  cin >> Q02 >> Q2;

  // Load evolution
  double Q0 = sqrt(Q02) - eps;
  double Q  = sqrt(Q2);
  APFEL::EvolveAPFEL(Q0,Q);

  cout << scientific << setprecision(5) << endl;
  // Tabulate PDFs for the LHA x values
  cout << "alpha_QCD(mu2F) = " << APFEL::AlphaQCD(Q) << endl;
  cout << endl;

  cout << "   x   " 
       << setw(11) << "    u-ubar    " 
       << setw(11) << "   d-dbar    " 
       << setw(11) << " 2(ubr+dbr)  " 
       << setw(11) << " c+cbar  " 
       << setw(11) << "     gluon   " << endl;

  cout << scientific;
  for (int i = 1; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
	 << setw(11) <<APFEL::xPDF(2,xlha[i]) - APFEL::xPDF(-2,xlha[i]) << "  "
	 << setw(11) <<APFEL::xPDF(1,xlha[i]) - APFEL::xPDF(-1,xlha[i]) << "  "
	 << setw(11) <<2*(APFEL::xPDF(-1,xlha[i]) + APFEL::xPDF(-2,xlha[i])) << "  "
	 << setw(11) <<APFEL::xPDF(4,xlha[i]) + APFEL::xPDF(-4,xlha[i]) << "  "
	 << setw(11) <<APFEL::xPDF(0,xlha[i])
	 << endl;
  cout << "      " << endl;

  return 0;
}
