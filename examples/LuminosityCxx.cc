#include <iostream>
#include <cmath>
#include <iomanip>
#include "APFEL/APFEL.h"
using namespace std;

int main()
{
  APFEL::SetPerturbativeOrder(1);
  APFEL::SetNumberOfGrids(3);
  APFEL::SetGridParameters(1,130,3,1e-9);
  APFEL::SetGridParameters(2,60,5,1e-1);
  APFEL::SetGridParameters(3,20,5,8e-1);

  // Initializes integrals on the grids
  APFEL::InitializeAPFEL();

  double S = pow(8e3,2);
  double Q0 = sqrt(2e0);
  double Qmin = 10;
  double Qmax = 6e3;
  int NMX  = 30;

  cout << "Luminosities at sqrt(S) = " << sqrt(S) << " GeV:\n" << endl;
  cout << "- gg luminosity " << endl;
      
  for (int i=1; i < 30; i++)
    {
      double Q = Qmin * pow(Qmax/Qmin,(double)(i-1)/(double)(NMX-1));
      APFEL::EvolveAPFEL(Q0,Q);
      cout << "MX = " << Q << " LUMI = " << APFEL::LUMI(0,0,S) << endl; 
    }

  cout << endl;
  return 0;
}
