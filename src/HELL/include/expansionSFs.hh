
using namespace std;


// RCvar=0   central
// RCvar=1   variation of gamma_plus (using linear as approximation for RC)
// RCvar=2   variation of gamma_qg   (using r -> as*b0)


namespace HELLx {

  double beta0(int nf);
  //double fmom(double N) { return 4*N/(N+1)/(N+1); }
  //
  // expansion of the LLp/NLL
  double Paux0(double x, int nf);
  double Paux0sq(double x, int nf);
  double Paux0cb(double x, int nf);
  double Paux1(double x, int nf, bool useLLp);
  double Paux0Paux1(double x, int nf, bool useLLp);
  double Paux2(double x, int nf, bool useLLp, int RCvar=0);

  // expansion of the LL
  double PLL1(double x, int nf);
  double PLL2(double x, int nf, int RCvar=0);

  // expansion of the NLL
  double PNLL2(double x, int nf);  // this is independent of var
  double PNLL3(double x, int nf, int RCvar=0);

  // momentum conservation functions
  double mcPgg1LL(double x, int nf);
  double mcPgg2LL(double x, int nf, int RCvar=0);
  double mcPgg2NLL(double x, int nf, bool useLLp);
  double mcPgg3NLL(double x, int nf, bool useLLp, int RCvar=0);

  // Pqg
  double Pqg2(double x, int nf, bool useLLp);
  double Pqg3(double x, int nf, bool useLLp, int RCvar=0);


};

