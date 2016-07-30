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

  // Evolve PDFs on the grids
  double Q02, Q2, eps = 1e-10;
  cout << "Enter initial and final scale in GeV^2" << endl;
  cin >> Q02 >> Q2;

  double Q0 = sqrt(Q02) - eps;
  double Q  = sqrt(Q2);

  APFEL::EvolveAPFEL(Q0,Q);

  double momsr = 0.0;
  for (int i = -6; i < 7; i++)
    momsr += APFEL::NPDF(i,2);
  
  momsr += APFEL::Ngamma(2);

  double uvsr = APFEL::NPDF(2,1) - APFEL::NPDF(-2,1);
  double dvsr = APFEL::NPDF(1,1) - APFEL::NPDF(-1,1);
  double svsr = APFEL::NPDF(3,1) - APFEL::NPDF(-3,1);

  cout << setprecision(15);
  cout << "Sum rules at Q = " << Q << " GeV:\n" << endl;
  cout << "- Momentum sum rule         = " << momsr << endl;
  cout << "- Up valence sum rule       = " << uvsr << endl;
  cout << "- Down valence sum rule     = " << dvsr << endl;
  cout << "- Strange valence sum rule  = " << svsr << endl;

  return 0;
}
