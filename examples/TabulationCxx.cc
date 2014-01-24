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
  //APFEL::SetPerturbativeOrder(2);
  //APFEL::SetPDFSet("MRST2004qed.LHgrid");
  //APFEL::EnableEvolutionOperator(true);
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
	 << setw(11) <<APFEL::xPDF(2,xlha[i]) - APFEL::xPDF(-2,xlha[i]) << "  "
	 << setw(11) <<APFEL::xPDF(1,xlha[i]) - APFEL::xPDF(-1,xlha[i]) << "  "
	 << setw(11) <<2*(APFEL::xPDF(-1,xlha[i]) + APFEL::xPDF(-2,xlha[i])) << "  "
	 << setw(11) <<APFEL::xPDF(4,xlha[i]) + APFEL::xPDF(-4,xlha[i]) << "  "
	 << setw(11) <<APFEL::xPDF(0,xlha[i]) << "  "
	 << setw(11) <<APFEL::xgamma(xlha[i]) << "  "
	 << endl;

  int n = 10;
  double *xext = new double[n+1];
  double ****M = new double***[14];

  for (int i = 0; i < 14; i++) {
    M[i] = new double**[14];
    for (int j = 0; j < 14; j++) {
      M[i][j] = new double*[n+1];
      for (int k = 0; k < n+1; k++)
	M[i][j][k] = new double[n+1];
    }    
  }
  
  APFEL::ExternalEvolutionOperator(Q0,Q,n,xext,M);

  for (int i = 0; i < 14; i++){
    for (int j = 0; j < 14; j++){
      for (int z = 0; z < n+1; z++){
	if (!M[i][j][z]) delete[] M[i][j][z];
      }
      if (!M[i][j]) delete[] M[i][j];
    }
    if (!M[i]) delete[] M[i];
  }

  delete[] xext;  
  delete[] M;
    
  return 0;
}
