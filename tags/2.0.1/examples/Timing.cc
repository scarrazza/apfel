#include <iostream>
#include <iomanip>
#include <cmath>
#include <sys/time.h>
#include <stdlib.h>
#include "APFEL/APFEL.h"
using namespace std;

class mytimer {

 private:
  timeval ti;
  timeval tf;

 public:
  mytimer() {};
  void start()   { gettimeofday(&ti, NULL); }
  void stop()    { gettimeofday(&tf, NULL); }
  double print() {
    cout << fixed;
    cout << setprecision(3);
    double tot = (tf.tv_sec - ti.tv_sec) * 1000.0;
    tot += (tf.tv_usec - ti.tv_usec) / 1000.0;    
    return tot;
  } 
};


int main(int argc, char** argv)
{
  int pt = 1;
  string theory = "QCD";
  if (argc > 2)
    {
      pt = atoi(argv[1]);
      theory.assign(argv[2]);
    }
  else
    {
      cout << "\nusage: Timing [Order 0,1,2] [Theory]\n" << endl;
      exit(-1);
    }

  // Activate some options
  APFEL::SetPerturbativeOrder(pt);
  APFEL::SetQLimits(1,1000);
  APFEL::SetTheory(theory);
  
  // Initializes integrals on the grids
  APFEL::InitializeAPFEL();

  // Load evolution

  mytimer t;
  double Q0 = 1.0;

  const int N = 100;
  for (int i = 0; i < N; i++)
    {
      double Q = exp(log(1.0)+i*(log(1000)-log(1.0))/N);
      t.start();
      APFEL::EvolveAPFEL(Q0,Q);
      t.stop();    
      cout << Q << "\t" << t.print()/1000 << endl;
    }

  return 0;
}
