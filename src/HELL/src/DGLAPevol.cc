/* ----------------------------------------------------

   It solves a (2-dimensional) matrix equation of the form

   d
   -- f(t,N) = Gamma(t,N) * f(t,N)
   dt

   The solution is given in terms of an evolutor U

   f(t,N) = U(t,t0,N) f(t0,N)

   with initial condition

   U(t0,t0,N) = 1   (identity in 2 dimension)


   !!!!! It is possible to use this solution also for
   an evolution in alpha_s.
   Indeed the equation is the same, provided we identify

   t -> alpha_s
   Gamma -> Gamma / beta(alpha_s)

   where beta(alpha_s) is the QCD beta-function, computed
   at the appropriate order.

   ----------------------------------------------------- */


#include "include/DGLAPevol.hh"


using namespace std;



DGLAPevol::DGLAPevol(sqmatrix<dcomplex> (*EvolutionMatrixGamma)(double t, dcomplex N)) {
  EvolMatrix = EvolutionMatrixGamma;
  id2 = sqmatrix<dcomplex>(1.,0.,0.,1.);
}

DGLAPevol::~DGLAPevol() {
  return;
}

void DGLAPevol::SetEvolMatrix(sqmatrix<dcomplex> (*EvolutionMatrixGamma)(double t, dcomplex N)) {
  EvolMatrix = EvolutionMatrixGamma;
}






sqmatrix<dcomplex> DGLAPevol::Uk(double t, double dt, dcomplex N) {
  //return id2 + EvolMatrix(t, N) * (dcomplex)dt;
  return id2 + EvolMatrix(t, N) * dt;
}
sqmatrix<dcomplex> DGLAPevol::Uk_exp(double t, double dt, dcomplex N) {
  sqmatrix<dcomplex> gam = EvolMatrix(t,N);
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
  return Mp * exp(dt*gp) + Mm * exp(dt*gm);
}

sqmatrix<dcomplex> DGLAPevol::evolU_linearized(double t, double t0, dcomplex N, int n) {
  sqmatrix<dcomplex> U = id2;
  double tk, dt = (t-t0)/(n+1.);
  for(int k=0; k<=n; k++) {
    tk = t0 + (k+0.5)*dt;
    U = Uk(tk, dt, N) * U;  // Left product (path-ordering)
  }
  return U;
}
sqmatrix<dcomplex> DGLAPevol::evolU(double t, double t0, dcomplex N, int n) {
  sqmatrix<dcomplex> U = id2;
  double tk, dt = (t-t0)/(n+1.);
  for(int k=0; k<=n; k++) {
    tk = t0 + (k+0.5)*dt;
    U = Uk_exp(tk, dt, N) * U;  // Left product (path-ordering)
  }
  return U;
}





