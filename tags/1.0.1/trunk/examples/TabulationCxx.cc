#include <iostream>
#include <iomanip>
#include <cmath>
#include "APFEL/APFEL.h"
using namespace std;

int main()
{
  double xlha[] = {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 
		1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
  
  // Activate some options
  APFEL::SetPerturbativeOrder(2);
  APFEL::SetPDFSet("MRST2004qed.LHgrid");
  
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
  cout << "alpha_QED(mu2F) = " << APFEL::AlphaQED(Q) << endl;
  cout << endl;

  cout << "   x   " 
       << setw(11) << "    u-ubar    " 
       << setw(11) << "   d-dbar    " 
       << setw(11) << "  2(ubr+dbr) " 
       << setw(11) << "  c+cbar " 
       << setw(11) << "  gluon " 
       << setw(11) << "    photon " << endl;

  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
	 << setw(11) << APFEL::xPDF(2,xlha[i]) - APFEL::xPDF(-2,xlha[i]) << "  "
	 << setw(11) <<APFEL::xPDF(1,xlha[i]) - APFEL::xPDF(-1,xlha[i]) << "  "
	 << setw(11) <<2*(APFEL::xPDF(-1,xlha[i]) + APFEL::xPDF(-2,xlha[i])) << "  "
	 << setw(11) <<APFEL::xPDF(4,xlha[i]) + APFEL::xPDF(-4,xlha[i]) << "  "
	 << setw(11) <<APFEL::xPDF(0,xlha[i]) << "  "
	 << setw(11) <<APFEL::xgamma(xlha[i]) << "  "
	 << endl;

  return 0;
}
