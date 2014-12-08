/* -----------------------------------------

   Test of deltaPres

   ----------------------------------------- */
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <sys/time.h> 

#include "HELL/include/hell.hh"


using namespace std;





int main (int argc, char* argv[]) {

  struct timeval t0, t1;
  gettimeofday(&t0,NULL);


  string prepath = "./data/";


  // working at NLL
  HELL::LogOrder order = HELL::NLL;

  // example with fixed nf
  int nf = 5;
  HELL::HELLnf sxD(nf, order, prepath);

  // alternatively for variable nf
  //HELL::HELL sxD(order, prepath);
  // and maybe you want to set mass thresholds different from default
  //init_as_thresholds  (asmc, asmb, asmt);



  HELL::Order fixed_order_to_be_matched_to = HELL::NNLO;
  double as = 0.12;
  double x = 0.01;
  sqmatrix<double> xdPNLL = sxD.DeltaP(as, x, fixed_order_to_be_matched_to);  // gives x*DeltaP(x,as)
  cout << endl << "Printing x*DeltaP_ij(x=" << x << ", as=" << as << ")"  << endl;
  cout << "gg: " << xdPNLL.gg() << endl
       << "gq: " << xdPNLL.gq() << endl
       << "qg: " << xdPNLL.qg() << endl
       << "qq: " << xdPNLL.qq() << endl;

  dcomplex N = 1.+I/2.;
  sqmatrix<dcomplex> dgammaNLL = sxD.DeltaGamma(as, N, fixed_order_to_be_matched_to); // in this notation N=0 is the rightmost pole.
  cout << endl << "Printing DeltaGamma_ij(N=" << N << ", as=" << as << ")" << endl;
  cout << "gg: " << dgammaNLL.gg() << endl
       << "gq: " << dgammaNLL.gq() << endl
       << "qg: " << dgammaNLL.qg() << endl
       << "qq: " << dgammaNLL.qq() << endl;



  // Finish time
  gettimeofday(&t1,NULL);
  double t=t1.tv_sec-t0.tv_sec+(t1.tv_usec-t0.tv_usec)*0.000001;
  cout << endl << "Total time: " << t << "s" << endl << endl;

  return 0;

}
