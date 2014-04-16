#include <iostream>
#include <iomanip>
#include <cmath>
#include "APFEL/APFEL.h"
using namespace std;

void print(double *);

int main()
{
  double xlha[] = {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 
		1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
  
  // Activate some options
  APFEL::SetQLimits(0.5,20000);
  APFEL::SetTheory("QavDP");      
  APFEL::SetPerturbativeOrder(1);
  
  // Initializes integrals on the grids
  APFEL::InitializeAPFEL();

  double Q02, Q2, eps = 1e-10;
  cout << "Enter initial and final scale in GeV^2" << endl;
  cin >> Q02 >> Q2;

  // Load evolution
  double Qi, Qf;
  double Q0 = sqrt(Q02) - eps;
  double Q  = sqrt(Q2);

  // QECD
  for (int n = 1; n <= 100; n+=99)
    {
      cout << "QECDP with " << n << " steps" << endl;
      cout << endl;

      APFEL::SetTheory("QECDP");
      APFEL::SetPDFSet("ToyLH");

      const double delta = log(Q/Q0)/ (double) n;
      Qi = Q0;
      for (int i = 1; i <= n; i++)
	{
	  Qf = Qi * exp(delta);
	  APFEL::EvolveAPFEL(Qi,Qf);
	  APFEL::SetPDFSet("apfel");
	  Qi = Qf;
	}

      print(xlha);
    }

  // QCED
  for (int n = 1; n <= 100; n+=99)
    {
      cout << "QCEDP with " << n << " steps" << endl;
      cout << endl;

      APFEL::SetTheory("QCEDP");
      APFEL::SetPDFSet("ToyLH");

      const double delta = log(Q/Q0)/ (double) n;
      Qi = Q0;
      for (int i = 1; i <= n; i++)
	{
	  Qf = Qi * exp(delta);
	  APFEL::EvolveAPFEL(Qi,Qf);
	  APFEL::SetPDFSet("apfel");
	  Qi = Qf;
	}

      print(xlha);
    }  

  // QavD
  for (int n = 1; n <= 100; n+=99)
    {
      cout << "QavDP with " << n << " steps" << endl;
      cout << endl;

      APFEL::SetTheory("QavDP");
      APFEL::SetPDFSet("ToyLH");

      const double delta = log(Q/Q0)/ (double) n;
      Qi = Q0;
      for (int i = 1; i <= n; i++)
	{
	  Qf = Qi * exp(delta);
	  APFEL::EvolveAPFEL(Qi,Qf);
	  APFEL::SetPDFSet("apfel");
	  Qi = Qf;
	}

      print(xlha);
    }

  // QUniD
  for (int n = 1; n <= 100; n+=99)
    {
      cout << "QUniD with " << n << " steps" << endl;
      cout << endl;

      APFEL::SetTheory("QUniD");
      APFEL::SetPDFSet("ToyLH");

      const double delta = log(Q/Q0)/ (double) n;
      Qi = Q0;
      for (int i = 1; i <= n; i++)
	{
	  Qf = Qi * exp(delta);
	  APFEL::EvolveAPFEL(Qi,Qf);
	  APFEL::SetPDFSet("apfel");
	  Qi = Qf;
	}

      print(xlha);
    }
  return 0;
}

void print(double *xlha)
{
  cout << scientific << setprecision(5);
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
  cout << endl;

}
