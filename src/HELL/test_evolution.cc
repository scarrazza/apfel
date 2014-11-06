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
double as(double t) {
  return as0/(1+as0*beta0(_nf)*t);
}


dcomplex psi(dcomplex Z){
  dcomplex SUB = 0. ;
  dcomplex ZZ = Z;
  if(abs(imag(ZZ))<10.) { // if too close to the real axis...
  label1:
    if(real(ZZ)<10.) { // ...use recurrence relation to push real(z) large enough
      SUB = SUB - 1./ ZZ;
      ZZ = ZZ + 1.;
      goto label1;
    }
  }
  dcomplex RZ = 1./ ZZ;
  dcomplex DZ = RZ * RZ;
  // SUB + asympt expansion (Abramowitz, Stengun, 6.3.18)
  return SUB + log(ZZ) - 0.5 * RZ - DZ/5040. * ( 420.+ DZ * ( - 42. + DZ * (20. - 21. * DZ) ) );
}
const double EulerGamma = 0.5772156649015328606065120900824;


// LO
dcomplex gamma0qq(dcomplex N) {
  return HELL::CF/2./M_PI * ( 1./N - 1./(N+1.) - 2.*psi(N+1.) - 2.*EulerGamma + 3./2. );
}
dcomplex gamma0qg(dcomplex N, double nf) {
  return nf/2./M_PI * ( 2.+N+N*N ) / N/(N+1.)/(N+2.);
}
dcomplex gamma0gq(dcomplex N) {
  return HELL::CF/2./M_PI * ( 2.+N+N*N ) / N/(N*N-1.);
}
dcomplex gamma0gg(dcomplex N, double nf) {
  return HELL::CA/M_PI * ( 11./12. + 1./(N-1.) - 1./N + 1./(N+1.) - 1./(N+2.) - psi(N+1.) - EulerGamma ) - nf/6./M_PI;
}
sqmatrix<dcomplex> gamma_LO(dcomplex N) {
  return sqmatrix<dcomplex>( gamma0gg(N,_nf), gamma0gq(N), gamma0qg(N,_nf), gamma0qq(N) );
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

// resummed evolution matrix
HELL::HELLnf *sxD;
sqmatrix<dcomplex> gammaRes(double t, dcomplex N) {
  return as(t) * gamma_LO(N) + sxD->DeltaGamma(as(t), N, HELL::LO);
}
sqmatrix<dcomplex> gammaRes_as(double as, dcomplex N) {
  return ( as*gamma_LO(N) + sxD->DeltaGamma(as,N,HELL::LO) ) / ( -beta0(_nf)*as*as );
}










int main() {

  //Start time counter
  struct timeval time0, time1;
  gettimeofday(&time0,NULL);


  _nf = 5;
  cout << endl << "Using  nf = " << _nf << endl << endl;


  double t = log(100.), t0=0.;
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
  PO.SetEvolMatrix(gammaRes);
  U = PO.evolU(t, t0, N, nsub);
  cout << "Resummed evolution at LO+NLL with  n = " << nsub << "  subdivisions:" << endl
       << U << endl;

  PO.SetEvolMatrix(gammaRes_as);
  U = PO.evolU(as(t), as(t0), N, nsub);
  cout << "Resummed evolution at LO+NLL with  n = " << nsub << "  subdivisions using as evolution:" << endl
       << U << endl;






  // Finish time
  gettimeofday(&time1,NULL);
  double deltatime=time1.tv_sec-time0.tv_sec+(time1.tv_usec-time0.tv_usec)*0.000001;
  cout << endl << "Total execution time: " << deltatime << "s" << endl;

  return 0;
}

