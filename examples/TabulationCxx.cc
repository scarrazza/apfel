#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "APFEL/APFEL.h"
using namespace std;

int main()
{
  double xlha[] = {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 
		1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
  
  // Activate some options
  //APFEL::SetFastEvolution(true);
  //APFEL::SetPerturbativeOrder(0);
  //APFEL::SetPDFSet("MRST2004qed");
  //APFEL::EnableEvolutionOperator(true);
  // Initializes integrals on the grids
  //APFEL::SetTheory("QCD");
  //APFEL::SetSmallxResummation(0, "NLL");
  //APFEL::SetPDFSet("NNPDF31_nnlo_as_0118");
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

  cout << "Standard evolution:" << endl;
  cout << "   x   " 
       << setw(11) << "    u-ubar    " 
       << setw(11) << "   d-dbar    " 
       << setw(11) << " 2(ubr+dbr)  " 
       << setw(11) << " c+cbar  " 
       << setw(11) << "     gluon   " 
       << setw(11) << "     photon  "
       << setw(11) << "    e^-+e^+  "
       << setw(11) << "   mu^-+mu^+ "
       << setw(11) <<"   tau^-+tau^+" << endl;

  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
	 << setw(11) << APFEL::xPDFj(2,xlha[i]) - APFEL::xPDFj(-2,xlha[i]) << "  "
	 << setw(11) << APFEL::xPDFj(1,xlha[i]) - APFEL::xPDFj(-1,xlha[i]) << "  "
	 << setw(11) << 2*(APFEL::xPDFj(-1,xlha[i]) + APFEL::xPDFj(-2,xlha[i])) << "  "
	 << setw(11) << APFEL::xPDFj(4,xlha[i]) + APFEL::xPDFj(-4,xlha[i]) << "  "
	 << setw(11) << APFEL::xPDFj(0,xlha[i]) << "  "
	 << setw(11) << APFEL::xgammaj(xlha[i]) << "  "
	 << setw(11) << APFEL::xLeptonj(1,xlha[i]) + APFEL::xLeptonj(-1,xlha[i]) << "  "
	 << setw(11) << APFEL::xLeptonj(2,xlha[i]) + APFEL::xLeptonj(-2,xlha[i]) << "  "
	 << setw(11) << APFEL::xLeptonj(3,xlha[i]) + APFEL::xLeptonj(-3,xlha[i]) << "  "
	 << endl;
  cout << "      " << endl;

  cout << "Standard evolution using the xPDFall function:" << endl;
  cout << "   x   " 
       << setw(11) << "    u-ubar    " 
       << setw(11) << "   d-dbar    " 
       << setw(11) << " 2(ubr+dbr)  " 
       << setw(11) << " c+cbar  " 
       << setw(11) << "     gluon   " << endl;

  cout << scientific;
  double xf[13];
  for (int i = 2; i < 11; i++) {
    APFEL::xPDFall(xlha[i],xf);
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
	 << setw(11) << xf[2+6] - xf[-2+6] << "  "
	 << setw(11) << xf[1+6] - xf[-1+6] << "  "
	 << setw(11) << 2*(xf[-1+6] + xf[-2+6]) << "  "
	 << setw(11) << xf[4+6] + xf[-4+6] << "  "
	 << setw(11) << xf[0+6] << "  "
	 << endl;
  }
  cout << "      " << endl;
  //
  // Cached PDFs
  //
  APFEL::CachePDFsAPFEL(Q0);

  cout << "Cached evolution:" << endl;
  cout << "   x   " 
       << setw(11) << "    u-ubar    " 
       << setw(11) << "   d-dbar    " 
       << setw(11) << " 2(ubr+dbr)  " 
       << setw(11) << " c+cbar  " 
       << setw(11) << "     gluon   " 
       << setw(11) << "     photon  "
       << setw(11) << "    e^-+e^+  "
       << setw(11) << "   mu^-+mu^+ "
       << setw(11) <<"   tau^-+tau^+" << endl;

  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
	 << setw(11) << APFEL::xPDFxQ(2,xlha[i],Q) - APFEL::xPDFxQ(-2,xlha[i],Q) << "  "
	 << setw(11) << APFEL::xPDFxQ(1,xlha[i],Q) - APFEL::xPDFxQ(-1,xlha[i],Q) << "  "
	 << setw(11) << 2*(APFEL::xPDFxQ(-1,xlha[i],Q) + APFEL::xPDFxQ(-2,xlha[i],Q)) << "  "
	 << setw(11) << APFEL::xPDFxQ(4,xlha[i],Q) + APFEL::xPDFxQ(-4,xlha[i],Q) << "  "
	 << setw(11) << APFEL::xPDFxQ(0,xlha[i],Q) << "  "
	 << setw(11) << APFEL::xPDFxQ(22,xlha[i],Q) << "  "
	 << setw(11) << APFEL::xPDFxQ(11,xlha[i],Q) + APFEL::xPDFxQ(-11,xlha[i],Q) << "  "
	 << setw(11) << APFEL::xPDFxQ(13,xlha[i],Q) + APFEL::xPDFxQ(-13,xlha[i],Q) << "  "
	 << setw(11) << APFEL::xPDFxQ(15,xlha[i],Q) + APFEL::xPDFxQ(-15,xlha[i],Q) << "  "
	 << endl;
  cout << "      " << endl;

  cout << "Cached evolution using the xPDFxQall function:" << endl;
  cout << "   x   " 
       << setw(11) << "    u-ubar    " 
       << setw(11) << "   d-dbar    " 
       << setw(11) << " 2(ubr+dbr)  " 
       << setw(11) << " c+cbar  " 
       << setw(11) << "     gluon   " << endl;

  cout << scientific;
  for (int i = 2; i < 11; i++) {
    APFEL::xPDFxQall(xlha[i],Q,xf);
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
	 << setw(11) << xf[2+6] - xf[-2+6] << "  "
	 << setw(11) << xf[1+6] - xf[-1+6] << "  "
	 << setw(11) << 2*(xf[-1+6] + xf[-2+6]) << "  "
	 << setw(11) << xf[4+6] + xf[-4+6] << "  "
	 << setw(11) << xf[0+6] << "  "
	 << endl;
  }
  cout << "      " << endl;

  return 0;
}
