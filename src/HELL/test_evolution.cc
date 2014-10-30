/* ----------------------------------------------------

   It solves a (2-dimensional) matrix equation of the form

   d
   -- f(t,N) = Gamma(t,N) * f(t,N)
   dt

   The solution is given in terms of an evolutor U

   f(t,N) = U(t,t0,N) f(t0,N)

   with initial condition

   U(t0,t0,N) = 1   (identity in 2 dimension)


   !!!!! It is possible to use this solution also for an evolution in alpha_s.
   Indeed the equation is the same, provided we identify

   t -> alpha_s
   Gamma -> Gamma / beta(alpha_s)

   where beta(alpha_s) is the QCD beta-function, computed at the appropriate order.

   ----------------------------------------------------- */

#include <iostream>
#include <sys/time.h>
#include "include/DGLAPevol.hh"
#include "include/math/matrix.hh"
#include "include/hell.hh"
//#include "../QCD.hh"
//#include "../anomalous_dimensions/gammaLO.hh"
//#include "../anomalous_dimensions/gammaNLO.hh"

using namespace std;


double _nf;
double as0 = 0.1;


double beta0(int nf) {
  return (11.*HELL::CA - 2.*nf)/12./M_PI;
}
double beta1(int nf){
  using namespace HELL;
  return (17.*CA*CA - (10.*CA*0.5+6.*CF*0.5)*nf)/(24.*M_PI*M_PI);
}
double as(double t) {
  return as0/(1.+as0*beta0(_nf)*t);
}




#include "include/gammaNLO.hh"
// LO
sqmatrix<dcomplex> gamma_LO(dcomplex N) {
  return sqmatrix<dcomplex>( gamma0gg(N,_nf), gamma0gq(N), gamma0qg(N,_nf), gamma0qq(N) );
}
// NLO
sqmatrix<dcomplex> gamma_NLO(dcomplex N) {
  gamma1sums g1s;
  sums(N, g1s);
  return sqmatrix<dcomplex>( gamma1SGgg(N,_nf,g1s), gamma1SGgq(N,_nf,g1s), gamma1SGqg(N,_nf,g1s), gamma1SGqq(N,_nf,g1s) );
}



sqmatrix<dcomplex> evolU_LO_exact(double t, double t0, dcomplex N) {
  double asint = -log(as(t)/as(t0))/beta0(_nf);  // integral of the LO as
  sqmatrix<dcomplex> gam = gamma_LO(N);
  dcomplex htr = gam.trace()/2.;
  dcomplex d = gam.det();
  dcomplex gp = htr + sqrt(htr*htr - d);
  dcomplex gm = htr - sqrt(htr*htr - d);
  dcomplex cp = gam.entry21() / (gp - gam.entry22());
  dcomplex cm = gam.entry21() / (gm - gam.entry22());
  sqmatrix<dcomplex> Mp(cm, -1., cp*cm, -cp);
  Mp /= (cm-cp);
  sqmatrix<dcomplex> Mm(-cp, 1., -cp*cm, cm);
  Mm /= (cm-cp);
  return Mp * exp(asint*gp) + Mm * exp(asint*gm);
}






// Fixed order evolution matrix
sqmatrix<dcomplex> gamma(double t, dcomplex N) {
  return as(t) * gamma_LO(N);
  // using evolution in as (then the variable t is in fact as)
  //return -( gamma_LO(N) ) / ( beta0(_nf)*t );
}
sqmatrix<dcomplex> gammaLO_as(double as, dcomplex N) {
  return ( as*gamma_LO(N) ) / ( -beta0(_nf)*as*as );
}
sqmatrix<dcomplex> gammaNLO_as(double as, dcomplex N) {
  return ( as*gamma_LO(N) + as*as*gamma_NLO(N) ) / ( -beta0(_nf)*as*as -beta1(_nf)*as*as*as);
}

// resummed evolution matrix
HELL::HELLnf *sxD;
sqmatrix<dcomplex> gammaResLO(double t, dcomplex N) {
  return as(t) * gamma_LO(N) + sxD->DeltaGamma(as(t), N, HELL::LO);
}
sqmatrix<dcomplex> gammaResLO_as(double as, dcomplex N) {
  return ( as*gamma_LO(N) + sxD->DeltaGamma(as,N,HELL::LO) ) / ( -beta0(_nf)*as*as );
}
sqmatrix<dcomplex> gammaResNLO_as(double as, dcomplex N) {
  return ( as*gamma_LO(N) + as*as*gamma_NLO(N) + sxD->DeltaGamma(as,N,HELL::NLO) ) / ( -beta0(_nf)*as*as -beta1(_nf)*as*as*as);
}










int main() {

  //Start time counter
  struct timeval time0, time1;
  gettimeofday(&time0,NULL);


  _nf = 5;
  cout << endl << "Using  nf = " << _nf << endl << endl;

  double Q = 10;
  double t = 2*log(Q), t0=0.;
  dcomplex N = I*4.+3.;
  int nsub = 100;

  DGLAPevol PO(gamma);

  sqmatrix<dcomplex> U, Ue, Udiff;
  U    = PO.evolU(t, t0, N, nsub);
  //U    = PO.evolU(as(t), as(t0), N, nsub);
  Ue   = evolU_LO_exact(t, t0, N);
  Udiff = Ue-U;


  cout << "Exact LO result:" << endl
       << Ue << endl;
  cout << "Using path-ordering with  n = " << nsub << "  subdivisions:" << endl
       << U << endl;
  cout << "Difference with exact result:" << endl
       << Udiff << endl << endl;




  // resummation part

  string prepath = "./data/";
  sxD = new HELL::HELLnf(_nf, HELL::NLL, prepath);
  PO.SetEvolMatrix(gammaResLO);
  //U = PO.evolU(t, t0, N, nsub);
  cout << "Resummed evolution at LO+NLL with  n = " << nsub << "  subdivisions:" << endl
       << U << endl;

  PO.SetEvolMatrix(gammaResLO_as);
  //U = PO.evolU(as(t), as(t0), N, nsub);
  cout << "Resummed evolution at LO+NLL with  n = " << nsub << "  subdivisions using as evolution:" << endl
       << U << endl;



  // test for Valerio
  cout << "--- TEST for Valerio ---" << endl;
  //
  // LO
  //PO.SetEvolMatrix(gammaLO_as);
  //PO.SetEvolMatrix(gammaNLO_as);
  PO.SetEvolMatrix(gammaResLO_as);
  //
  as0 =  0.25538582434016766;
  double asc =  0.25000000000000000;
  double asb =  0.18637330935258775;
  double as1 =  0.10829923787758310;
  _nf=5;
  as0 = 0.24;
  as1 = as(log(10000./2.));
  cout << as1 << endl;
  //
  N = 2;
  /*
  _nf = 3;
  delete sxD;
  sxD = new HELL::HELLnf(_nf, HELL::NLL, prepath);
  U = PO.evolU(asc, as0, N, nsub);
  _nf = 4;
  delete sxD;
  sxD = new HELL::HELLnf(_nf, HELL::NLL, prepath);
  U = PO.evolU(asb, asc, N, nsub) * U;
  _nf = 5;
  delete sxD;
  sxD = new HELL::HELLnf(_nf, HELL::NLL, prepath);
  U = PO.evolU(as1, asb, N, nsub) * U;
  */
  delete sxD;
  sxD = new HELL::HELLnf(_nf, HELL::NLL, prepath);
  U = PO.evolU(as1, as0, N, nsub);
  //U = PO.evolU(t, t0, N, nsub);
  //
  //double g0 = 0.364857610;
  //double s0 = 0.635141627;
  //vec2<dcomplex> f0(g0, s0);
  //vec2<dcomplex> f0(0.36479821649757477, 0.63510921112928409);
  vec2<dcomplex> f0(0.364857610, 0.635141627);
  //vec2<dcomplex> f0(0.003729506,  0.027428291);
  vec2<dcomplex> f1 = U*f0;
  //vec2<dcomplex> f2 = Ue*f0;
  //
  cout << "initial condition:" << endl << f0 << endl;
  cout << "evolved:" << endl << f1 << endl;
  //cout << "evolved LO exact:" << endl << f2 << endl;




  // NLO
  PO.SetEvolMatrix(gammaNLO_as);
  PO.SetEvolMatrix(gammaResNLO_as);
  //
  as0 =  0.25617614682386308;
  asc =  0.25000000000000000;
  asb =  0.18157141297389037;
  as1 =  0.10431617707311804;
  //
  _nf = 3;
  delete sxD;
  sxD = new HELL::HELLnf(_nf, HELL::NLL, prepath);
  U = PO.evolU(as0, asc, N, nsub);
  _nf = 4;
  delete sxD;
  sxD = new HELL::HELLnf(_nf, HELL::NLL, prepath);
  U = PO.evolU(asc, asb, N, nsub) * U;
  _nf = 5;
  delete sxD;
  sxD = new HELL::HELLnf(_nf, HELL::NLL, prepath);
  U = PO.evolU(asb, as1, N, nsub) * U;
  //
  f1 = U*f0;
  //
  cout << f1 << endl;





  // Finish time
  gettimeofday(&time1,NULL);
  double deltatime=time1.tv_sec-time0.tv_sec+(time1.tv_usec-time0.tv_usec)*0.000001;
  cout << endl << "Total execution time: " << deltatime << "s" << endl;

  return 0;
}

