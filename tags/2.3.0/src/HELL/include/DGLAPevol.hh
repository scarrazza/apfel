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


#ifndef _PATH_ORDERING_
#define _PATH_ORDERING_

#include "./math/matrix.hh"


using namespace std;




class DGLAPevol {

private:
  sqmatrix<dcomplex> (*EvolMatrix)(double t, dcomplex N);
  sqmatrix<dcomplex> id2;
  sqmatrix<dcomplex> Uk(double t, double dt, dcomplex N);
  sqmatrix<dcomplex> Uk_exp(double t, double dt, dcomplex N);

public:
  DGLAPevol(sqmatrix<dcomplex> (*EvolutionMatrixGamma)(double t, dcomplex N));
  ~DGLAPevol();
  void SetEvolMatrix(sqmatrix<dcomplex> (*EvolutionMatrixGamma)(double t, dcomplex N));
  sqmatrix<dcomplex> evolU(double t, double t0, dcomplex N, int n);
  sqmatrix<dcomplex> evolU_linearized(double t, double t0, dcomplex N, int n);

};




#endif
