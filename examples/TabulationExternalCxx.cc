#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "APFEL/APFEL.h"
using namespace std;

extern "C" void externalsetapfel_(double const& x, double const& Q, double* xf);

int main()
{
  double xlha[] = {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 
		1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
  
  // Activate some options
  APFEL::SetPDFSet("external");
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
       << endl;

  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
	 << setw(11) << APFEL::xPDFj(2,xlha[i]) - APFEL::xPDFj(-2,xlha[i]) << "  "
	 << setw(11) << APFEL::xPDFj(1,xlha[i]) - APFEL::xPDFj(-1,xlha[i]) << "  "
	 << setw(11) << 2*(APFEL::xPDFj(-1,xlha[i]) + APFEL::xPDFj(-2,xlha[i])) << "  "
	 << setw(11) << APFEL::xPDFj(4,xlha[i]) + APFEL::xPDFj(-4,xlha[i]) << "  "
	 << setw(11) << APFEL::xPDFj(0,xlha[i]) << "  "
	 << setw(11) << APFEL::xgammaj(xlha[i]) << "  "
	 << endl;
  cout << "      " << endl;

  return 0;
}

void externalsetapfel_(double const& x, double const& Q, double* xf)
{
  const double N_uv = 5.107200;
  const double auv  = 0.8;
  const double buv  = 3;
  const double N_dv = 3.064320;
  const double adv  = 0.8;
  const double bdv  = 4;
  const double N_g  = 1.7;
  const double ag   = -0.1;
  const double bg   = 5;
  const double N_db = 0.1939875;
  const double adb  = -0.1;
  const double bdb  = 6;
  const double fs   = 0.2;

  // User defined PDFs
  const double xuv   = N_uv * pow(x, auv) * pow(1 - x, buv);
  const double xdv   = N_dv * pow(x, adv) * pow(1 - x, bdv);
  const double xg    = N_g  * pow(x, ag)  * pow(1 - x, bg);
  const double xdbar = N_db * pow(x, adb) * pow(1 - x, bdb);
  const double xubar = xdbar * ( 1 - x );
  const double xs    = fs * ( xdbar + xubar );
  const double xsbar = xs;

  // Initialize PDFs to zero
  for (int i = -6; i <= 7; i++)
    xf[i+6] = 0;

  if (x > 1)
    return;

  xf[3+6]  = xs;
  xf[2+6]  = xuv + xubar;
  xf[1+6]  = xdv + xdbar;
  xf[0+6]  = xg;
  xf[-1+6] = xdbar;
  xf[-2+6] = xubar;
  xf[-3+6] = xsbar;
}
