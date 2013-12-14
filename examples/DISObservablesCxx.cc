#include <iostream>
#include <string>
#include <iomanip>
#include <cmath>
#include "APFEL/APFEL.h"
using namespace std;

int main()
{

  double q20 = 2.0;
  double q2  = 10;
  double y   = 0.5;
  double pol = 0.0;
  double xin = 1e-5;
  double xfi = 1e-1;
  int np  = 5;
  double q0  = sqrt(q20);
  double q   = sqrt(q2);
  string proj   = "ELECTRON";
  string target = "PROTON";
  int pto    = 2;
  string scheme = "FONLL";
  const int irep   = 0;
  string pdfset = "toyLH_NNLO.LHgrid";
  double F2[5],F3[5],FL[5],sigma[5];

  string proc = "EM";
  cout << "   x   " 
       << setw(11) << "  Q2    " 
       << setw(11) << "   y    " 
       << setw(11) << "   F2l  " 
       << setw(11) << "   F2c  " << endl;
  cout << scientific;


  double x = xin;
  for (int ip=1; ip <=np ;ip++)
    {
      APFEL::DIS_xsec(x,q0,q,y,pol,proc,scheme,pto,pdfset,irep,target,
                      proj,F2,F3,FL,sigma);

      cout << x << "\t" << q2 << "\t" << y << "\t" 
	   << F2[0] << "\t" << F2[1] << endl;

      x = x * exp( 1.0 / ( np - 1.0 ) * log( xfi / xin ) );
    }

  proc = "CC";
cout << "\n   x   " 
       << setw(11) << "  Q2    " 
       << setw(11) << "   y    " 
       << setw(11) << "   F2c  " 
       << setw(11) << "   FLc  " 
       << setw(11) << "   F3c  "
       << setw(11) << " sigmac " << endl;
  cout << scientific;


  x = xin;
  for (int ip=1; ip <=np ;ip++)
    {
      APFEL::DIS_xsec(x,q0,q,y,pol,proc,scheme,pto,pdfset,irep,target,
                      proj,F2,F3,FL,sigma);

      cout << x << "\t" << q2 << "\t" << y << "\t" 
	   << F2[1] << "\t" << FL[1] << "\t" << F3[1] << "\t" 
	   << sigma[1] << endl;

      x = x * exp( 1.0 / ( np - 1.0 ) * log( xfi / xin ) );
    }

  

  return 0;
}
