/* -----------------------------------------

   _  _ ____ _    _    
   |__| |___ |    |    
   |  | |___ |___ |___ x
   
   HELLx: High-Energy Large Logarithms - fast x-space version

   Author: Marco Bonvini

   Computes
   Delta P_res = P_res - P_FixedOrder
   Delta C_res = C_res - C_FixedOrder
   as an interpolation from pre-prepared input files

   ----------------------------------------- */

#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>

#include "HELL/include/hell-x.hh"
#include "HELL/include/math/special_functions.hh"


using namespace std;



namespace HELLx {

  bool useLLp = false;
  void SetLLpMode(bool use) {
    useLLp = use;
  }

  int _damping = 2, _dampingsqrt = 4;
  double beta0(int nf) { return (11*CA-2*nf)/12./M_PI; }
  double fmom(double N) { return 4*N/(N+1)/(N+1); }
  //
  const double c0 = CA/M_PI*4*log(2.);
  const double k0 = CA/M_PI*28.*ZETA3;
  // from DL-LO kernel
  inline double c1LO (int nf) { return -15.00496429 - 0.04503163717*nf;}
  inline double k1LO (int nf) { return -507.744719  - 1.080759292*nf;}
  // from DL-NLO kernel
  inline double c1NLO(int nf) { return -15.04043249 - 0.2076545137*nf; }
  inline double k1NLO(int nf) { return -574.496771  - 0.3724353405*nf; }
  inline double m1NLO(int nf) { return 0.3099881585 - 0.01878716112*nf; }
  //
  inline double a11()       { return CA/M_PI; }
  inline double a10(int nf) { return -(11.*CA + 2.*nf*(1.-2.*CF/CA))/12./M_PI; }
  inline double a21(int nf) { return nf*(26*CF-23*CA)/36/M_PI/M_PI; }
  //
  unsigned int factorials[] = {1,1,2,6,24,120,720};
  unsigned int factorial(unsigned int k) {
    if(k<7) return factorials[k];
    double res = 1;
    for(int j=1; j<=k; j++) res *= j;
    return res;
  }
  double binomial(unsigned int k, unsigned int j) {
    return factorial(k)/(1.*factorial(k-j)*factorial(j));
  }
  double Pole(double x, int k, int j) { // inverse Mellin of 1/N^k/(1+N)^j
    if(j==0 && k>0) return pow(-log(x),k-1)/factorial(k-1)/x;
    if(k==0 && j>0) return pow(-log(x),j-1)/factorial(j-1);
    if(k==1 && j==1) return 1/x-1;
    if(k==2 && j==1) return Pole(x,2,0) - Pole(x,1,0) + Pole(x,0,1);
    if(k==1 && j==2) return Pole(x,1,0) - Pole(x,0,2) - Pole(x,0,1);
    cout << "HELLx warning: this inverse mellin is not implemented: 1/N^"<<k<<"/(1+N)^"<<j << endl;
    return 0;
  }
  double Poly(double x, int k, int j) { // inverse Mellin of 1/N^k/(1+N)^j (psi_1(1+N)-Zeta2), up to O(x^0)
    double res = 0;
    if     (k==2 && j==0) res = -2*ZETA3/x + 2 - log(x);
    else if(k==3 && j==0) res = 3*ZETA4/x +2*ZETA3*log(x)/x -3 +log(x);
    else if(k==0 && j==2) res = -2*ZETA3 - pow(log(x),3)/6;
    else if(k==0 && j==3) res = 3*ZETA4 + 2*ZETA3*log(x) + pow(log(x),4)/24;
    else cout << "HELLx warning: this inverse mellin is not implemented: (psi_1(1+N)-Zeta2)/N^"<<k<<"/(1+N)^"<<j << endl;
    return res;
  }
  double exactPoly(double x, int k, int j) { // inverse Mellin of 1/N^k/(1+N)^j (psi_1(1+N)-Zeta2)
    double res = 0;
    if     (k==2 && j==0) res = (2*(Li3(x)-ZETA3) - Li2(x)*log(x))/x;
    else if(k==3 && j==0) res = (3*(ZETA4-Li4(x))+(Li3(x)+2*ZETA3)*log(x))/x;
    else if(k==0 && j==2) res = 2*(Li3(x)-ZETA3) - Li2(x)*log(x) - pow(log(x),3)/6;
    else if(k==0 && j==3) res = 3*(ZETA4-Li4(x))+(Li3(x)+2*ZETA3)*log(x) + pow(log(x),4)/24;
    else cout << "HELLx warning: this inverse mellin is not implemented: (psi_1(1+N)-Zeta2)/N^"<<k<<"/(1+N)^"<<j << endl;
    return res;
  }
  double d1(int k, int j) {
    double res = 0;
    for(int i=0; i<=k; i++) {
      res += binomial(k,i) * pow(-1.,i) /(2.+i) /binomial(4+2*i+j,j);
    }
    return res;
  }
  // expansion of the LLp
  double PLLp0(double x, int nf) {
    double a1 = a11(), a0=a10(nf);
    double res = a1/x + a0;
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k) * pow(1-sqrt(x),j);
  }
  double PLLp0sq(double x, int nf) {
    double a1 = a11(), a0=a10(nf);
    double res = a1*a1*Pole(x,2,0) + a0*a0*Pole(x,0,2) + 2*a1*a0*Pole(x,1,1);
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k) * pow(1-sqrt(x),j);
  }
  double PLLp0cb(double x, int nf) {
    double a1 = a11(), a0=a10(nf);
    double res = a1*a1*a1*Pole(x,3,0) + a0*a0*a0*Pole(x,0,3) + 3*a1*a1*a0*Pole(x,2,1) + 3*a1*a0*a0*Pole(x,1,2);
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k) * pow(1-sqrt(x),j);
  }
  double PLLp1(double x, int nf) {
    double res = 0;
    if(useLLp) {
      res = beta0(nf)*(3*k0/32.-c0) * (1/x-4*(1+log(x)));
    } else {
      double a1 = a21(nf);
      res = a1 * (1/x-2);
    }
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k) * pow(1-sqrt(x),j);
  }
  double PLLp0PLLp1(double x, int nf) {
    double a1 = a11(), a0=a10(nf);
    double res = 0;
    if(useLLp) {
      res = beta0(nf)*(3*k0/32.-c0) * ( a1*(Pole(x,2,0)-4*Pole(x,0,2)) + a0*(Pole(x,1,0)+4*Pole(x,0,3)-4*Pole(x,0,2)-Pole(x,0,1)) );
    } else {
      res = a21(nf) * ( a1*Pole(x,2,0) + (a0-2*a1)/x -2*a0*Pole(x,0,2) + 2*a1-a0 );
    }
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k) * pow(1-sqrt(x),j);
  }
  double PLLp2(double x, int nf, int var=0) {
    double res = 0;
    double a1 = a11(), a0=a10(nf);
    double chi02 = a1*a0;
    double b0 = beta0(nf);
    if(useLLp) {
      res = a1*chi02 * ( Pole(x,2,0) - 3*Pole(x,1,0) + 2*Pole(x,0,2) + 3*Pole(x,0,1) )
	+ 2*(a0 + a1)*chi02 * (2*Pole(x,0,4)-Pole(x,0,3));
      res += CA/M_PI* ( a1*a1*  exactPoly(x,3,0) 
			+a1*a0* exactPoly(x,2,0)
			+4*pow(a0 + a1,2)*              exactPoly(x,0,3)
			-2*(a0*a0 + 4*a0*a1 + 3*a1*a1)* exactPoly(x,0,2)
			);
      // add RC
      int T = (var==1 ? 1 : 2);
      double dc1 = c1NLO(nf)-c1LO(nf);
      double dk1 = k1NLO(nf)-k1LO(nf);
      double  m1 = m1NLO(nf);
      double rho2 = b0*(48*c0*k0-3*k0*k0-256*c0*c0)/256. + (16*c0*(dk1+6*k0*m1)+k0*(16*dc1-3*dk1-15*k0*m1))/512;
      double rho1 = b0*(3*T*k1NLO(nf)/32 -T*c1NLO(nf) + (b0+6*m1)*k0/16);
      //double  rho2 = -0.189252486317084 + 0.198339561094134*nf;
      //double  rho1 = -63.1502484908168 + 3.83820436998940*nf - 0.000661611070737904*nf*nf;
      //if(var) rho1 = -29.1702672639587 + 1.62760436906141*nf + 0.00850246161412463*nf*nf;
      res += rho1*Pole(x,1,0) + rho2*Pole(x,2,0) - 4*(1+log(x))*(rho1+rho2);
    } else {
      double A = ( CA*CA*(-74./27. +11./12.*ZETA2 +5./2.*ZETA3) +nf*((4*CA+7*CF)/27. +(CA-2*CF)*ZETA2/6.) )/M_PI/M_PI;
      double chi11 = a21(nf);
      res = (A*a1 + a1*chi02 + a1*chi11) * Pole(x,2,0)
	+ (-4*A*a1 - 3*a1*chi02 - 3*a1*chi11) * Pole(x,1,0)
	+ 4*(a0 + a1)*(2*A + chi02 + chi11) * Pole(x,0,4)
	- 2*(a0 + a1)*(4*A + chi02 + chi11) * Pole(x,0,3)
	+ 2*(A*a0 + 3*A*a1 + a1*chi02 + a1*chi11) * Pole(x,0,2)
	+ (4*A*a1 + 3*a1*chi02 + 3*a1*chi11) * Pole(x,0,1)
	;
      res += CA/M_PI* ( a1*a1*  exactPoly(x,3,0) 
			+a1*a0* exactPoly(x,2,0)
			+4*pow(a0 + a1,2)*              exactPoly(x,0,3)
			-2*(a0*a0 + 4*a0*a1 + 3*a1*a1)* exactPoly(x,0,2)
			);
      // add RC
      res += b0*b0*k0/16* (1/x-4*(1+log(x)));
    }
    // result
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k)*pow(1-sqrt(x),j);
  }
  // expansion of the LL
  double PLL1(double x, int nf) { return PLLp1(x,nf); }
  double PLL2(double x, int nf, int var=0) {
    double res = 0;
    double a1 = a11(), a0=a10(nf);
    double chi02 = a1*a0;
    res = a1*chi02 * ( Pole(x,2,0) - 3*Pole(x,1,0) + 2*Pole(x,0,2) + 3*Pole(x,0,1) )
      + 2*(a0 + a1)*chi02 * (2*Pole(x,0,4)-Pole(x,0,3));
    res += CA/M_PI* ( a1*a1*  exactPoly(x,3,0) 
		      +a1*a0* exactPoly(x,2,0)
		      +4*pow(a0 + a1,2)*              exactPoly(x,0,3)
		      -2*(a0*a0 + 4*a0*a1 + 3*a1*a1)* exactPoly(x,0,2)
		      );
    // add RC
    int T = (var==1 ? 1 : 2);
    double b0 = beta0(nf);
    double rho2 = b0*(48*c0*k0-3*k0*k0-256*c0*c0)/256.;
    double rho1 = b0*(3*T*k1LO(nf)/32 -T*c1LO(nf) + b0*k0/16);
    //double  rho2 = -2.76616449350881 + 0.167646332939928*nf;
    //double  rho1 = -55.5269225000529 + 3.17343574292416*nf + 0.0116261992352771*nf*nf;
    //if(var) rho1 = -26.9938523096279 + 1.49343193929256*nf + 0.00863994604701791*nf*nf;
    res += rho1*Pole(x,1,0) + rho2*Pole(x,2,0) - 4*(1+log(x))*(rho1+rho2);
    // result
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k)*pow(1-sqrt(x),j);
  }
  // expansion of the NLL
  double PNLL2(double x, int nf, int var=0) {  // this is independent of var
    double res = 0;
    double b0 = beta0(nf);
    double A = ( CA*CA*(-74./27. +11./12.*ZETA2 +5./2.*ZETA3) +nf*((4*CA+7*CF)/27. +(CA-2*CF)*ZETA2/6.) )/M_PI/M_PI;
    double a1 = a11(), a0=a10(nf);
    double chi11 = a21(nf);
    double chi02 = a1*a0;
    res = (A*a1 + a1*chi02 + a1*chi11) * Pole(x,2,0)
      + (-4*A*a1 - 3*a1*chi02 - 3*a1*chi11) * Pole(x,1,0)
      + 4*(a0 + a1)*(2*A + chi02 + chi11) * Pole(x,0,4)
      - 2*(a0 + a1)*(4*A + chi02 + chi11) * Pole(x,0,3)
      + 2*(A*a0 + 3*A*a1 + a1*chi02 + a1*chi11) * Pole(x,0,2)
      + (4*A*a1 + 3*a1*chi02 + 3*a1*chi11) * Pole(x,0,1)
      ;
    res += CA/M_PI* ( a1*a1*  Poly(x,3,0) 
		      +a1*a0* Poly(x,2,0)
		      +4*pow(a0 + a1,2)*              Poly(x,0,3)
		      -2*(a0*a0 + 4*a0*a1 + 3*a1*a1)* Poly(x,0,2)
		      );
    // add RC
    res += b0*b0*k0/16* (1/x-4*(1+log(x)));
    // result
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k)*pow(1-sqrt(x),j);
  }
  double PNLL3(double x, int nf, int var=0) {
    return 0; // to be written
    double res = 0;
    double b0 = beta0(nf);
    double A = ( CA*CA*(-74./27. +11./12.*ZETA2 +5./2.*ZETA3) +nf*((4*CA+7*CF)/27. +(CA-2*CF)*ZETA2/6.) )/M_PI/M_PI;
    double a1 = a11(), a0=a10(nf);
    double chi11 = a21(nf);
    double chi02 = a1*a0;
    res = (A*a1 + a1*chi02 + a1*chi11) * Pole(x,2,0)
      + (-4*A*a1 - 3*a1*chi02 - 3*a1*chi11) * Pole(x,1,0)
      + 4*(a0 + a1)*(2*A + chi02 + chi11) * Pole(x,0,4)
      - 2*(a0 + a1)*(4*A + chi02 + chi11) * Pole(x,0,3)
      + 2*(A*a0 + 3*A*a1 + a1*chi02 + a1*chi11) * Pole(x,0,2)
      + (4*A*a1 + 3*a1*chi02 + 3*a1*chi11) * Pole(x,0,1)
      ;
    res += CA/M_PI* ( a1*a1*  Poly(x,3,0) 
		      +a1*a0* Poly(x,2,0)
		      +4*pow(a0 + a1,2)*              Poly(x,0,3)
		      -2*(a0*a0 + 4*a0*a1 + 3*a1*a1)* Poly(x,0,2)
		      );
    // add RC
    res += b0*b0*k0/16* (1/x-4*(1+log(x)));
    // result
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k)*pow(1-sqrt(x),j);
  }
  // anomalous dimensions, needed for momcons
  double gammagg1LL(double N, int nf) {
    return beta0(nf)* (3*k0/32.-c0) * (1./N - fmom(N));
  }
  double gammagg2LL(double N, int nf, int var) {
    double res = 0;
    double a1 = a11(), a0=a10(nf);
    double chi02 = a1*a0;
    double gamma0  = a1/N + a0 - 2*(a1+a0)*N/(N+1);
    double gamma0p = -a1/N/N - 2*(a1+a0)/(N+1)/(N+1);
    res  = -gamma0p * ( chi02/pow(N+1.,2) - chi02/4*fmom(N) + a1*gamma0*(dpsi(N+1)-ZETA2) );
    // add RC
    int T = (var==1 ? 1 : 2);
    double b0 = beta0(nf);
    double rho2 = b0*(48*c0*k0-3*k0*k0-256*c0*c0)/256.;
    double rho1 = b0*(3*T*k1LO(nf)/32 -T*c1LO(nf) + b0*k0/16);
    res += rho2/N/N + rho1/N - (rho2+rho1)*fmom(N);
    // result
    return res;
  }
  double gammagg2NLL(double N, int nf, int var) { // this is independent of var
    double res = 0;
    double a1 = a11(), a0=a10(nf);
    double A = ( CA*CA*(-74./27. +11./12.*ZETA2 +5./2.*ZETA3) +nf*((4*CA+7*CF)/27. +(CA-2*CF)*ZETA2/6.) )/M_PI/M_PI;
    //double chi11 = a10(nf)*a11();
    //double chi02 = a21(nf);
    double chi11 = a21(nf);
    double chi02 = a1*a0;
    double gamma0p = -a1/N/N - 2*(a1+a0)/(N+1)/(N+1);
    res  = -gamma0p * ( A + chi11/(N+1.) + chi02/pow(N+1.,2) - (A+chi11/2.+chi02/4.)*fmom(N) );
    res += CA/M_PI* ( a1*a1*  ( 3*ZETA4/N -2*ZETA3/N/N -3/(N+1) -1/(N+1)/(N+1) )
		      +a1*a0* (-2*ZETA3/N +2/(N+1) +1/(N+1)/(N+1) )
		      +4*pow(a0 + a1,2)*              ( 3*ZETA4/(N+1) -2*ZETA3/(N+1)/(N+1) +1/pow(N+1,5) )
		      -2*(a0*a0 + 4*a0*a1 + 3*a1*a1)* ( -2*ZETA3/(N+1) +1/pow(N+1,4) )
		      );
    // add RC
    double b0 = beta0(nf);
    res += b0*b0*k0/16.* (1/N - fmom(N));
    // add Pqg
    double gamma0LL = a1/N + a0/(N+1);
    double gamma1LL = (useLLp ? b0*(3*k0/32.-c0) * (1/N - fmom(N)) : a21(nf) * (1/N-2/(N+1)) );
    double h1 = 5./3., h2 = 14./9.;
    res += (1-CF/CA) * (h2*gamma0LL*(gamma0LL-b0) + h1*gamma1LL) * nf/3./M_PI;
    // result
    return res;
  }
  double gammagg3NLL(double N, int nf, int var) {
    return 0;  // to be written
    /*
    double res = 0;
    double a1 = a11(), a0=a10(nf);
    // add RC
    double b0 = beta0(nf);
    res += 0;
    // add Pqg
    double gamma0LL = a1/N + a0/(N+1);
    double gamma1LL = (useLLp ? b0*(3*k0/32.-c0) * (1/N - fmom(N)) : a21(nf) * (1/N-2/(N+1)) );
    double gamma2LL = 0;//(useLLp ? b0*(3*k0/32.-c0) * (1/N - fmom(N)) : a21(nf) * (1/N-2/(N+1)) );
    double h1 = 5./3., h2 = 14./9., h3 = 2*(41./81.+ZETA3);
    //res += (1-CF/CA) * (h2*gamma0LL*(gamma0LL-b0) + h1*gamma1LL) * nf/3./M_PI;
    // result
    return res;
    */
  }
  // momentum conservation functions
  double mcPgg1LL(double x, int nf) {
    double momcons = 0;
    int k=_damping, j=_dampingsqrt;
    for(int n=0; n<=k; n++) {
      for(int m=0; m<=j; m++) {
	momcons += binomial(k,n) * binomial(j,m) * pow(-1.,n+m) * gammagg1LL(1.+n+m/2.,nf);
      }
    }
    momcons /= d1(k,j);
    return momcons * pow(1-x,k)*pow(1-sqrt(x),j);
  }
  double mcPgg2LL(double x, int nf, int var=0) {
    double momcons = 0;
    int k=_damping, j=_dampingsqrt;
    if(_damping==2 && _dampingsqrt==4) {
      double coeffs[9] = {1, -4, 4, 4, -10, 4, 4, -4, 1};
      for(int n=0; n<=8; n++) {
	momcons += coeffs[n] * gammagg2LL(1.+n/2.,nf,var);
      }
    } else {
      for(int n=0; n<=k; n++) {
	for(int m=0; m<=j; m++) {
	  momcons += binomial(k,n) * binomial(j,m) * pow(-1.,n+m) * gammagg2LL(1.+n+m/2.,nf,var);
	}
      }
    }
    momcons /= d1(k,j);
    return momcons * pow(1-x,k)*pow(1-sqrt(x),j);
  }
  double mcPgg2NLL(double x, int nf, int var=0) {
    double momcons = 0;
    int k=_damping, j=_dampingsqrt;
    if(_damping==2 && _dampingsqrt==4) {
      double coeffs[9] = {1, -4, 4, 4, -10, 4, 4, -4, 1};
      for(int n=0; n<=8; n++) {
	momcons += coeffs[n] * gammagg2NLL(1.+n/2.,nf,var);
      }
    } else {
      for(int n=0; n<=k; n++) {
	for(int m=0; m<=j; m++) {
	  momcons += binomial(k,n) * binomial(j,m) * pow(-1.,n+m) * gammagg2NLL(1.+n+m/2.,nf,var);
	}
      }
    }
    momcons /= d1(k,j);
    return momcons * pow(1-x,k)*pow(1-sqrt(x),j);
  }
  double mcPgg3NLL(double x, int nf, int var=0) {
    double momcons = 0;
    int k=_damping, j=_dampingsqrt;
    if(_damping==2 && _dampingsqrt==4) {
      double coeffs[9] = {1, -4, 4, 4, -10, 4, 4, -4, 1};
      for(int n=0; n<=8; n++) {
	momcons += coeffs[n] * gammagg3NLL(1.+n/2.,nf,var);
      }
    } else {
      for(int n=0; n<=k; n++) {
	for(int m=0; m<=j; m++) {
	  momcons += binomial(k,n) * binomial(j,m) * pow(-1.,n+m) * gammagg3NLL(1.+n+m/2.,nf,var);
	}
      }
    }
    momcons /= d1(k,j);
    return momcons * pow(1-x,k)*pow(1-sqrt(x),j);
  }
  // Pqg
  double Pqg2(double x, int nf) {
    double h1 = 5./3., h2 = 14./9.;
    double b0 = beta0(nf);
    double P0sq = PLLp0sq(x,nf);
    double P0   = PLLp0  (x,nf);
    double P1   = PLLp1  (x,nf);
    double res = h2*(P0sq-b0*P0) + h1*P1;
    return res*nf/3./M_PI;
  }
  double Pqg3(double x, int nf, int RCvar) {
    double h1 = 5./3., h2 = 14./9., h3 = 2*(41./81.+ZETA3);
    double b0 = beta0(nf);
    double P0cb = PLLp0cb(x,nf);
    double P0sq = PLLp0sq(x,nf);
    double P0   = PLLp0  (x,nf);
    double P1   = PLLp1  (x,nf);
    double P2   = PLLp2  (x,nf); //,RCvar);  // here we don;t have to vary it, as in Pqg we always use the same gammaLLp
    double P0P1 = PLLp0PLLp1(x,nf);
    double T = (RCvar==1 ? 1 : 2);
    double res = h3*(P0cb-3*b0*P0sq+2*b0*b0*P0) + h2*(2*P0P1-T*b0*P1) + h1*P2;
    return res*nf/3./M_PI;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////







  double Qofalphas(double as, double as0, double Q0, int nf) {
    double b0 = (33 - 2*nf)/12./M_PI;
    double b1 = (153-19*nf)/24./M_PI/M_PI/b0;
    return Q0 * exp(1./2./b0 * (1./as - 1./as0 - b1*log((b1 + 1./as)/(b1 + 1./as0))));
  }
  const double mZ = 91.1876;
  const double asZ[4] = {0.10466, 0.11227, 0.118, 0.118};
  double Qofalphas(double as, int nf) {
    return Qofalphas(as, asZ[nf-3], mZ, nf);
  }





  double minterpolate(double mQ, double *mQvals, double *F, int Nmass, double x, double as, int nf) {
    if(mQ>mQvals[Nmass-1]) {
      cout << "\033[0;31m" << "HELLx: Warning! Extrapolating out of interpolation range: m/Q=" << mQ << " > " << mQvals[Nmass-1] << " for as=" << as << ", nf=" << nf << "\033[0m" << endl;
    }
    if(mQ<mQvals[0]) {
      cout << "\033[0;31m" << "HELLx: Warning! Extrapolating out of interpolation range: m/Q=" << mQ << " < " << mQvals[0]       << " for as=" << as << ", nf=" << nf << "\033[0m" << endl;
    }
    int m0 = -1;
    for(int m=0; m<Nmass; m++) {
      if(mQ > mQvals[m]) m0++;
    }
    if(m0 == -1)      m0 = 0;
    if(m0 == Nmass-1) m0 = Nmass-2;
    //if(m0 == -1)      return __builtin_nan("");
    //if(m0 == Nmass-1) return __builtin_nan("");
    //if(m0 == -1)      return F[0];
    //if(m0 == Nmass-1) return F[Nmass-1];
    double resLin = F[m0] + (F[m0+1]-F[m0])/(mQvals[m0+1]-mQvals[m0])*(mQ-mQvals[m0]);
    double resLog = F[m0]*exp(log(F[m0+1]/F[m0])/log(mQvals[m0+1]/mQvals[m0])*log(mQ/mQvals[m0]));
    if(F[m0]<=0 || F[m0+1]<=0) return resLin;
    double w = log(1/x);
    return (resLin*w+resLog)/(1+w);
  }





  xTable::xTable(string filename) {
    infile = new ifstream(filename.c_str());
    //ifstream infile(filename.c_str());
    if(!infile->good()) {
      cout << "\033[0;31m" << "HELLx: Error reading table " << filename << "\033[0m" << endl;
      exit(0);
    }
    *infile >> Np1 >> Np2 >> x_min >> x_mid >> x_max;
    xx    = new double[Np1+Np2];
    for(int i=0; i<Np1+Np2; i++) {
      if(i<Np1) xx[i] = x_min * exp(i/(Np1-1.)*log(x_mid/x_min));
      else      xx[i] = x_mid + (i-Np1+1)*(x_max-x_mid)/(Np2-0.);
    }
    //xx[Np1+Np2] = 1.;
    //infile.close();
  }
  void xTableP::Init() {
    xdPgg = new double[Np1+Np2];
    xdPqg = new double[Np1+Np2];
    for(int i=0; i<Np1+Np2; i++) {
      *infile >> xdPgg[i];
      if(isNLL) *infile >> xdPqg[i];
    }
    infile->close();
  }
  void xTableC::Init() {
    xdC2g = new double[Np1+Np2];
    xdCLg = new double[Np1+Np2];
    for(int i=0; i<Np1+Np2; i++) {
      *infile >> xdC2g[i] >> xdCLg[i];
    }
    infile->close();
  }
  void xTableCm::Init() {
    *infile >> Nmass >> Q;
    mQvals  = new double [Nmass];
    xdKhg   = new double*[Nmass];
    xdC2g   = new double*[Nmass];
    xdCLg   = new double*[Nmass];
    xdC2axg = new double*[Nmass];
    xdCLaxg = new double*[Nmass];
    xdC2CCg = new double*[Nmass];
    xdCLCCg = new double*[Nmass];
    xdC3CCg = new double*[Nmass];
    for(int m=0; m<Nmass; m++) {
      xdKhg  [m] = new double[Np1+Np2];
      xdC2g  [m] = new double[Np1+Np2];
      xdCLg  [m] = new double[Np1+Np2];
      xdC2axg[m] = new double[Np1+Np2];
      xdCLaxg[m] = new double[Np1+Np2];
      xdC2CCg[m] = new double[Np1+Np2];
      xdCLCCg[m] = new double[Np1+Np2];
      xdC3CCg[m] = new double[Np1+Np2];
      *infile >> mQvals[m];
      for(int i=0; i<Np1+Np2; i++) {
	*infile >> xdKhg[m][i] >> xdC2g[m][i] >> xdCLg[m][i] >> xdC2axg[m][i] >> xdCLaxg[m][i] >> xdC2CCg[m][i] >> xdCLCCg[m][i] >> xdC3CCg[m][i];
      }
    }
    infile->close();
  }
  void xTableCggH::Init() {
    xdCggH    = new double[Np1+Np2];
    xdCggHaux = new double[Np1+Np2];
    *infile >> c10 >> c20 >> c11 >> c30 >> c21;
    for(int i=0; i<Np1+Np2; i++) {
      *infile >> xdCggH[i] >> xdCggHaux[i];
    }
    infile->close();
  }
  double xTable::interpolate(double x) {
    if(x>1 || x<0) {
      cout << "\033[0;31m" << "HELLx: Error! Requesting resummed splitting function for unphysical value of x=" << x << " outside the physical range 0<x<=1" << "\033[0m" << endl;
      exit(45);
    }
    if(x>x_max) {
      cout << "\033[0;31m" << "HELLx: Warning! Extrapolating out of interpolation range: x=" << x << " > x_max=" << x_max << "\033[0m" << endl;
      x = x_max;
    }
    if(x<x_min) {
      cout << "\033[0;31m" << "HELLx: Warning! Extrapolating out of interpolation range: x=" << x << " < x_min=" << x_min << "\033[0m" << endl;
      x = x_min;
    }
    double ii;
    if(x<x_mid) ii = (Np1-1.)*log(x/x_min)/log(x_mid/x_min);
    else        ii = Np1-1.+Np2*(x-x_mid)/(x_max-x_mid);
    return ii;
  }
  void xTableP::eval(double x, double &dPgg, double &dPqg) {
    double ii = interpolate(x);
    int i = int(ii);
    if(i<0) {
      cout << "\033[0;31m" << "HELLx: Error! This should never happen" << "\033[0m" << endl;
      abort();
    }
    dPgg = xdPgg[i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdPgg[i+1]-xdPgg[i]));
    dPqg = 0;
    if(isNLL)
      dPqg = xdPqg[i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdPqg[i+1]-xdPqg[i]));
    return;
  }
  void xTableC::eval(double x, double &dC2g, double &dCLg) {
    double ii = interpolate(x);
    int i = int(ii);
    if(i<0) {
      cout << "\033[0;31m" << "HELLx: Error! This should never happen" << "\033[0m" << endl;
      abort();
    }
    dC2g = xdC2g[i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdC2g[i+1]-xdC2g[i]));
    dCLg = xdCLg[i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdCLg[i+1]-xdCLg[i]));
    return;
  }
  void xTableCm::eval(double x, double mQ, double &dKhg, double &dC2g, double &dCLg, double &dC2axg, double &dCLaxg, double &dC2CCg, double &dCLCCg, double &dC3CCg, double as, int nf) {
    double ii = interpolate(x);
    int i = int(ii);
    if(i<0 || ii<i) {
      cout << "\033[0;31m" << "HELLx: Error! This should never happen" << "\033[0m" << endl;
      abort();
    }
    double *mdKhg, *mdC2g, *mdCLg, *mdC2axg, *mdCLaxg, *mdC2CCg, *mdCLCCg, *mdC3CCg;
    mdKhg   = new double[Nmass];
    mdC2g   = new double[Nmass];
    mdCLg   = new double[Nmass];
    mdC2axg = new double[Nmass];
    mdCLaxg = new double[Nmass];
    mdC2CCg = new double[Nmass];
    mdCLCCg = new double[Nmass];
    mdC3CCg = new double[Nmass];
    for(int m=0; m<Nmass; m++) {
      mdKhg  [m] = xdKhg  [m][i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdKhg  [m][i+1]-xdKhg  [m][i]));
      mdC2g  [m] = xdC2g  [m][i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdC2g  [m][i+1]-xdC2g  [m][i]));
      mdCLg  [m] = xdCLg  [m][i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdCLg  [m][i+1]-xdCLg  [m][i]));
      mdC2axg[m] = xdC2axg[m][i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdC2axg[m][i+1]-xdC2axg[m][i]));
      mdCLaxg[m] = xdCLaxg[m][i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdCLaxg[m][i+1]-xdCLaxg[m][i]));
      mdC2CCg[m] = xdC2CCg[m][i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdC2CCg[m][i+1]-xdC2CCg[m][i]));
      mdCLCCg[m] = xdCLCCg[m][i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdCLCCg[m][i+1]-xdCLCCg[m][i]));
      mdC3CCg[m] = xdC3CCg[m][i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdC3CCg[m][i+1]-xdC3CCg[m][i]));
    }
    dKhg   = minterpolate(mQ, mQvals, mdKhg,   Nmass, x, as, nf);
    dC2g   = minterpolate(mQ, mQvals, mdC2g,   Nmass, x, as, nf);
    dCLg   = minterpolate(mQ, mQvals, mdCLg,   Nmass, x, as, nf);
    dC2axg = minterpolate(mQ, mQvals, mdC2axg, Nmass, x, as, nf);
    dCLaxg = minterpolate(mQ, mQvals, mdCLaxg, Nmass, x, as, nf);
    dC2CCg = minterpolate(mQ, mQvals, mdC2CCg, Nmass, x, as, nf);
    dCLCCg = minterpolate(mQ, mQvals, mdCLCCg, Nmass, x, as, nf);
    dC3CCg = minterpolate(mQ, mQvals, mdC3CCg, Nmass, x, as, nf);
    delete[] mdKhg;
    delete[] mdC2g;
    delete[] mdCLg;
    delete[] mdC2axg;
    delete[] mdCLaxg;
    delete[] mdC2CCg;
    delete[] mdCLCCg;
    delete[] mdC3CCg;
    return;
  }
  void xTableCggH::eval(double x, double &dCggH, double &dCggHaux) {
    double ii = interpolate(x);
    int i = int(ii);
    if(i<0) {
      cout << "\033[0;31m" << "HELLx: Error! This should never happen" << "\033[0m" << endl;
      abort();
    }
    dCggH    = xdCggH   [i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdCggH   [i+1]-xdCggH   [i]));
    dCggHaux = xdCggHaux[i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdCggHaux[i+1]-xdCggHaux[i]));
    return;
  }








  // **************************
  // ***   nf fixed part    ***
  // **************************

  void HELLxnf::Init(string datapath) {
    datapath_ = datapath;
    //string pyexec = "python "+datapath+"/gen_info.py "+datapath;
    //system(pyexec.c_str());
    string sord = (_order == 1 ? "nlo" : "lo");
    ostringstream filename;
    filename << datapath << "/" << sord << "_nf" << _nf << ".info";
    ifstream info(filename.str().c_str());
    if(!info.good()) {
      cout << "\033[0;31m" << "HELLx: Error reading info file" << "\033[0m" << endl;
      cout << "Do you have the tables properly installed?" << endl
	   << "The latest set of tables can be downloaded from the webpage https://www.ge.infn.it/~bonvini/hell/ and must be placed in HELLx/data" << endl;
      cout << "If you are using HELLx through APFEL, place the tables in <APFELdir>/src/HELL/data and run 'make install' again" << endl;
      exit(0);
    }
    double astmp;
    while(info.good()) {
      info >> astmp;
      _alphas.push_back(astmp);
      // hack!!!
      if(astmp<=0.2) _alphasHgg.push_back(astmp);
    }
    _alphas.pop_back();
    info.close();
  }

  template<class S>
  void deleteTable(map<int,S*> &T) {
    typename map<int,S*>::iterator itxT;  // typename tells the compiler that what follows is a typename, and not a field
    for (itxT=T.begin(); itxT!=T.end(); ++itxT) {
      delete T[itxT->first];
    }
  }
  HELLxnf::~HELLxnf() {
    deleteTable(xT[0]);
    deleteTable(xT[1]);
    deleteTable(xTC[0]);
    deleteTable(xTC[1]);
    deleteTable(xTCm[0]);
    deleteTable(xTCm[1]);
    deleteTable(xTCggH[0]);
    deleteTable(xTCggH[1]);
  }


  int HELLxnf::GetOrder() {
    return _order;
  }
  int HELLxnf::GetNf() {
    return _nf;
  }

  void HELLxnf::GetAvailableAlphas(vector<double> &as) {
    as.resize(_alphas.size());
    as = _alphas;
  }


  int HELLxnf::alphas_interpolation(double as, vector<double> vas, double &factor) {
    int k = 0;
    if(as < vas[0] || as > vas[vas.size()-1]) {
      cout << "\033[0;31m" << "HELLx: ERROR: alpha_s=" << as << " out of interpolation range [" << vas[0] << ", " << vas[vas.size()-1] << "] for nf=" << _nf << "\033[0m" << endl;
      exit(22);
    }
    for(unsigned int i=1; i<vas.size(); i++) {
      //if(as <= vas[i]) break;
      //k++;
      if(as < vas[i]) {
	k = i-1;
	break;
      }
    }
    // the following is needed for cubic interpolation
    if(k==0) k++;
    if(k==vas.size()-1) k -= 2;
    if(k==vas.size()-2) k -= 1;
    //cout << as << "  " << k << endl;
    //
    //cout << alphas[k] << " < " << as << " < " << alphas[k+1] << endl;
    factor = (as-vas[k]) / (vas[k+1]-vas[k]);
    return k;
  }

  string sas(double as) {
    ostringstream os;
    if     (as<0.01) os << "000" << int(1000*as);
    else if(as<0.1 ) os << "00"  << int(1000*as);
    else if(as<1.  ) os << "0"   << int(1000*as);
    else os << int(1000*as);
    return os.str();
  }

  template<class S>
  void HELLxnf::ReadTable(int k, map<int,S*> &T, string basename) {
    string sord = (_order == 1 ? "nlo" : "lo");
    ostringstream filename;
    typename map<int,S*>::iterator itxT;  // typename tells the compiler that what follows is a typename, and not a field
    //for(int i=k; i<k+2; i++) { // for linear interpolation
    for(int i=k-1; i<k+3; i++) { // for cubic interpolation
      itxT = T.find(i);
      if (itxT == T.end()) {
	filename.str("");  filename.clear();
	filename << datapath_ << "/" << basename << "_" << sord << "_nf" << _nf << "_alphas" << sas(_alphas[i]) << ".table";
	T[i] = new S(filename.str(), _order);
      }
    }
  }
  double HELLxnf::alphas_cubicinterpolate(double a, double k, vector<double> vas, double *y) {
    double as[4] = { vas[k-1], vas[k], vas[k+1], vas[k+2] };
    //return (a-as[2])*y[1]/(as[1]-as[2]) + (a-as[1])*y[2]/(as[2]-as[1]);
    return (a-as[1])*(a-as[2])*(a-as[3])*y[0]/(as[0]-as[1])/(as[0]-as[2])/(as[0]-as[3])
      +    (a-as[0])*(a-as[2])*(a-as[3])*y[1]/(as[1]-as[0])/(as[1]-as[2])/(as[1]-as[3])
      +    (a-as[0])*(a-as[1])*(a-as[3])*y[2]/(as[2]-as[0])/(as[2]-as[1])/(as[2]-as[3])
      +    (a-as[0])*(a-as[1])*(a-as[2])*y[3]/(as[3]-as[0])/(as[3]-as[1])/(as[3]-as[2]);
  }

  sqmatrix<double> HELLxnf::DeltaP(double as, double x, Order matched_to_fixed_order) {
    if(_order==1 && (matched_to_fixed_order < NLO || matched_to_fixed_order > N3LO)) {
      cout << "Error: trying to match NLL resummation in splitting functions to a fixed order different from NLO or NNLO" << endl;
      exit(45);
    }
    if(_order==0 && (matched_to_fixed_order < LO || matched_to_fixed_order > NLO)) {
      cout << "Error: trying to match LL resummation in splitting functions to a fixed order different from LO or NLO" << endl;
      exit(45);
    }
    double factor;
    int k = alphas_interpolation(as, _alphas, factor);
    ReadTable(k, xT[_RCvar], "xtable");
    //double xdPgg[2], xdPqg[2];
    //xT[_RCvar][k  ]->eval(x,xdPgg[0],xdPqg[0]);
    //xT[_RCvar][k+1]->eval(x,xdPgg[1],xdPqg[1]);
    //double dGgg = (xdPgg[0]+factor*(xdPgg[1]-xdPgg[0]))/x;
    //double dGqg = (xdPqg[0]+factor*(xdPqg[1]-xdPqg[0]))/x;
    double xdPgg[4], xdPqg[4];
    xT[_RCvar][k-1]->eval(x,xdPgg[0],xdPqg[0]);
    xT[_RCvar][k  ]->eval(x,xdPgg[1],xdPqg[1]);
    xT[_RCvar][k+1]->eval(x,xdPgg[2],xdPqg[2]);
    xT[_RCvar][k+2]->eval(x,xdPgg[3],xdPqg[3]);
    double dGgg = alphas_cubicinterpolate(as, k, _alphas, xdPgg)/x;
    double dGqg = alphas_cubicinterpolate(as, k, _alphas, xdPqg)/x;
    if(_order==0 && matched_to_fixed_order == LO) {
      dGgg += as*as * ( PLL1(x,_nf) - mcPgg1LL(x,_nf) );
    }
    if(_order==0 && matched_to_fixed_order == NNLO) {
      dGgg -= as*as*as * ( PLL2(x,_nf,_RCvar) - mcPgg2LL(x,_nf,_RCvar) );
    }
    if(_order==1 && matched_to_fixed_order == NLO) {
      dGgg += as*as*as * ( PNLL2(x,_nf) - CF/CA*Pqg2(x,_nf) - mcPgg2NLL(x,_nf) );
      dGqg += as*as*as *   Pqg2(x,_nf);
      //dGqg += as*as*as*as *   Pqg3(x,_nf,_RCvar);
    }
    if(_order==1 && matched_to_fixed_order == NNLO) {
      //dGqg += as*as*as*as *   Pqg3(x,_nf,_RCvar);
    }
    if(_order==1 && matched_to_fixed_order == N3LO) {
      dGgg -= as*as*as*as * ( PNLL3(x,_nf,_RCvar) - CF/CA*Pqg3(x,_nf,_RCvar) - mcPgg3NLL(x,_nf,_RCvar) );
      dGqg -= as*as*as*as *   Pqg3(x,_nf,_RCvar);
    }
    return sqmatrix<double>(dGgg, CF/CA*dGgg, dGqg, CF/CA*dGqg);
  }
  double HELLxnf::DeltaC(double as, double x, Order matched_to_fixed_order, string id) {
    if(_order==0) return 0;
    if(matched_to_fixed_order < NLO || matched_to_fixed_order > N3LO) {
      cout << "Error: trying to match resummed massless coefficient functions to a fixed order different from NLO, NNLO or N3LO" << endl;
      exit(46);
    }
    double factor;
    int k = alphas_interpolation(as, _alphas, factor);
    ReadTable(k, xTC[_RCvar], "xtableC");
    //double xdC2g[2], xdCLg[2];
    //xTC[_RCvar][k  ]->eval(x,xdC2g[0],xdCLg[0]);
    //xTC[_RCvar][k+1]->eval(x,xdC2g[1],xdCLg[1]);
    double xdC2g[4], xdCLg[4];
    xTC[_RCvar][k-1]->eval(x,xdC2g[0],xdCLg[0]);
    xTC[_RCvar][k  ]->eval(x,xdC2g[1],xdCLg[1]);
    xTC[_RCvar][k+1]->eval(x,xdC2g[2],xdCLg[2]);
    xTC[_RCvar][k+2]->eval(x,xdC2g[3],xdCLg[3]);
    map<string,double> dCg;
    //dCg["F2"] = (xdC2g[0]+factor*(xdC2g[1]-xdC2g[0]))/x;
    //dCg["FL"] = (xdCLg[0]+factor*(xdCLg[1]-xdCLg[0]))/x;
    dCg["F2"] = alphas_cubicinterpolate(as, k, _alphas, xdC2g)/x;
    dCg["FL"] = alphas_cubicinterpolate(as, k, _alphas, xdCLg)/x;
    if(matched_to_fixed_order == NLO) {
      double P0 = PLLp0(x,_nf);
      double h12 = 43./18.-ZETA2;
      double h1L = -1./3;
      dCg["F2"] += as*as * P0 * h12*_nf/3./M_PI;
      dCg["FL"] += as*as * P0 * h1L*_nf/3./M_PI;
    }
    if(matched_to_fixed_order == N3LO) {
      double P0sq = PLLp0sq(x,_nf);
      double P0   = PLLp0  (x,_nf);
      double P1   = PLLp1  (x,_nf);
      double h12 = 43./18.-ZETA2;
      double h1L = -1./3;
      double h22 = 3.2982926655872 /2.;  // Q0MSbar
      //double h22 = 9.7092628157717 /2.;  // MSbar
      double h2L = 2.132843710929551495;
      double b0 = beta0(_nf);
      dCg["F2"] -= as*as*as * (h22*(P0sq-b0*P0) + h12*P1) *_nf/3./M_PI;
      dCg["FL"] -= as*as*as * (h2L*(P0sq-b0*P0) + h1L*P1) *_nf/3./M_PI;
    }
    return dCg[id];
  }
  double HELLxnf::DeltaCm(double as, double x, Order matched_to_fixed_order, string id, double mQ) {
    if(_order==0) return 0;
    if(matched_to_fixed_order < NLO || matched_to_fixed_order > NNLO) {
      cout << "HELLx error: trying to match resummed massive coefficient functions or matching conditions to a fixed order different from NLO or NNLO" << endl;
      exit(47);
    }
    double factor;
    int k = alphas_interpolation(as, _alphas, factor);
    //cout << _alphas[k] << "  " << _alphas[k+1] << "  " << factor << endl;
    ReadTable(k, xTCm[_RCvar], "xtableCm");
    // converts m/Q-interpolation into m-interpolation
    //double mQ0 = mQ * Qofalphas(as,_nf)/Qofalphas(_alphas[k  ],_nf);
    //double mQ1 = mQ * Qofalphas(as,_nf)/Qofalphas(_alphas[k+1],_nf);
    double mQ0 = mQ * Qofalphas(as,_nf)/Qofalphas(_alphas[k-1],_nf);
    double mQ1 = mQ * Qofalphas(as,_nf)/Qofalphas(_alphas[k  ],_nf);
    double mQ2 = mQ * Qofalphas(as,_nf)/Qofalphas(_alphas[k+1],_nf);
    double mQ3 = mQ * Qofalphas(as,_nf)/Qofalphas(_alphas[k+2],_nf);
    //double xdKhg[2], xdC2g[2], xdCLg[2];
    //xTCm[_RCvar][k  ]->eval(x, mQ0, xdKhg[0],xdC2g[0],xdCLg[0],as,_nf);
    //xTCm[_RCvar][k+1]->eval(x, mQ1, xdKhg[1],xdC2g[1],xdCLg[1],as,_nf);
    double xdKhg[4], xdC2g[4], xdCLg[4], xdC2axg[4], xdCLaxg[4], xdC2CCg[4], xdCLCCg[4], xdC3CCg[4];
    xTCm[_RCvar][k-1]->eval(x, mQ0, xdKhg[0],xdC2g[0],xdCLg[0],xdC2axg[0],xdCLaxg[0],xdC2CCg[0],xdCLCCg[0],xdC3CCg[0],as,_nf);
    xTCm[_RCvar][k  ]->eval(x, mQ1, xdKhg[1],xdC2g[1],xdCLg[1],xdC2axg[1],xdCLaxg[1],xdC2CCg[1],xdCLCCg[1],xdC3CCg[1],as,_nf);
    xTCm[_RCvar][k+1]->eval(x, mQ2, xdKhg[2],xdC2g[2],xdCLg[2],xdC2axg[2],xdCLaxg[2],xdC2CCg[2],xdCLCCg[2],xdC3CCg[2],as,_nf);
    xTCm[_RCvar][k+2]->eval(x, mQ3, xdKhg[3],xdC2g[3],xdCLg[3],xdC2axg[3],xdCLaxg[3],xdC2CCg[3],xdCLCCg[3],xdC3CCg[3],as,_nf);
    //cout << setw(14) << mQ
    //	 << setw(14) << xdKhg[0] << setw(14) << xdC2g[0] << setw(14) << xdCLg[0]
    //	 << setw(14) << xdKhg[1] << setw(14) << xdC2g[1] << setw(14) << xdCLg[1] << endl;
    map<string,double> dCg;
    //dCg["Khg"] = (xdKhg[0]+factor*(xdKhg[1]-xdKhg[0]))/x;
    //dCg["F2m"] = (xdC2g[0]+factor*(xdC2g[1]-xdC2g[0]))/x;
    //dCg["FLm"] = (xdCLg[0]+factor*(xdCLg[1]-xdCLg[0]))/x;
    dCg["Khg"  ] = alphas_cubicinterpolate(as, k, _alphas, xdKhg)/x;
    dCg["F2m"  ] = alphas_cubicinterpolate(as, k, _alphas, xdC2g)/x;
    dCg["FLm"  ] = alphas_cubicinterpolate(as, k, _alphas, xdCLg)/x;
    dCg["F2axm"] = alphas_cubicinterpolate(as, k, _alphas, xdC2axg)/x;
    dCg["FLaxm"] = alphas_cubicinterpolate(as, k, _alphas, xdCLaxg)/x;
    dCg["F2CCm"] = alphas_cubicinterpolate(as, k, _alphas, xdC2CCg)/x;
    dCg["FLCCm"] = alphas_cubicinterpolate(as, k, _alphas, xdCLCCg)/x;
    dCg["F3CCm"] = alphas_cubicinterpolate(as, k, _alphas, xdC3CCg)/x;
    if(matched_to_fixed_order == NLO) {
      double P0 = PLLp0(x,_nf);
      double mQ2 = mQ*mQ;
      double lmQ2 = log(mQ2);
      double l1mQ2 = log(1+mQ2);
      double sq = sqrt(1+4*mQ2);
      double h1K   = -(28.+30.*lmQ2+9.*lmQ2*lmQ2)/18.;
      double h12   = (5+3*lmQ2)/6 - (1-mQ2)*HPLmp(-1./sq)/sq + ArcCsch(2*sqrt(mQ2))*(13-10*mQ2+6*(1-mQ2)*lmQ2)/3/sq;
      double h1L   = (-1+12*mQ2+3*(1+6*mQ2)*lmQ2)/3/(1+4*mQ2) + 4*mQ2*(1+3*mQ2)*HPLmp(-1./sq)/(1+4*mQ2)/sq + ArcCsch(2*sqrt(mQ2))*(6+8*mQ2*(1-6*mQ2)-24*mQ2*(1+3*mQ2)*lmQ2)/3/(1+4*mQ2)/sq;
      double h12ax = 12*mQ2*(8*log(4*mQ2) - 16*log(1+sq) - 3*HPLmp(-1./sq) + ArcCsch(2*sqrt(mQ2))*(6*lmQ2+28))/6/sq;
      double h1Lax = 2*mQ2*(2*ArcCsch(2*sqrt(mQ2))*(28+82*mQ2+6*(7*mQ2+2)*lmQ2) - sq*(8+3*log(mQ2)) - 6*mQ2*log((sq-1)/(sq+1)) - 6*(2+7*mQ2)*HPLmp(-1./sq))/3/(1+4*mQ2)/sq;
      double h12CC = (58-30*lmQ2-9*lmQ2*lmQ2+78*l1mQ2+18*l1mQ2*l1mQ2-36*Li2(1/(1+mQ2)))/36;
      double h1LCC = (4*(-3+16*mQ2+9*l1mQ2)+3*mQ2*(-28*mQ2*ArcCoth(1+2*mQ2)+lmQ2*(-2+(3+6*mQ2)*lmQ2)+30*l1mQ2-6*mQ2*l1mQ2*l1mQ2)+36*mQ2*mQ2*Li2(1/(1+mQ2)))/36/(1+mQ2);
      double h13CC = (28+30*lmQ2+9*(1+2*mQ2)*lmQ2*lmQ2+60*mQ2*(lmQ2-l1mQ2)-18*mQ2*l1mQ2*l1mQ2+36*mQ2*Li2(1/(1+mQ2)))/36./(1+mQ2);
      dCg["Khg"  ] += as*as * P0 * h1K  /3./M_PI;
      dCg["F2m"  ] += as*as * P0 * h12  /3./M_PI;
      dCg["FLm"  ] += as*as * P0 * h1L  /3./M_PI;
      dCg["F2axm"] += as*as * P0 * h12ax/3./M_PI;
      dCg["FLaxm"] += as*as * P0 * h1Lax/3./M_PI;
      dCg["F2CCm"] += as*as * P0 * h12CC/3./M_PI;
      dCg["FLCCm"] += as*as * P0 * h1LCC/3./M_PI;
      dCg["F3CCm"] += as*as * P0 * h13CC/3./M_PI;
    }
    return dCg[id];
  }
  double HELLxnf::DeltaCggH(double as, double x, double muFrat, Order matched_to_fixed_order, bool aux) {
    if(_order==0) return 0;
    if(matched_to_fixed_order < NLO || matched_to_fixed_order > N3LO) {
      cout << "Error: trying to match resummed massless coefficient functions to a fixed order different from NLO or NNLO" << endl;
      exit(46);
    }
    double factor;
    //int k = alphas_interpolation(as, _alphas, factor);
    int k = alphas_interpolation(as, _alphasHgg, factor);
    ostringstream name;
    name << "xtableChiggs" << (_RCvar==1 ? "var" : "") << 1000*muFrat;
    ReadTable(k, xTCggH[_RCvar], name.str());
    double xdCggH[4], xdCggHaux[4];
    xTCggH[_RCvar][k-1]->eval(x,xdCggH[0],xdCggHaux[0]);
    xTCggH[_RCvar][k  ]->eval(x,xdCggH[1],xdCggHaux[1]);
    xTCggH[_RCvar][k+1]->eval(x,xdCggH[2],xdCggHaux[2]);
    xTCggH[_RCvar][k+2]->eval(x,xdCggH[3],xdCggHaux[3]);
    double res = 0;
    if(aux) res = alphas_cubicinterpolate(as, k, _alphasHgg, xdCggHaux)/x;
    else    res = alphas_cubicinterpolate(as, k, _alphasHgg, xdCggH)/x;
    //if(matched_to_fixed_order >= NNLO) {
    if(matched_to_fixed_order < NNLO) {
      double P0sq = PLLp0sq(x,_nf);
      double P0   = PLLp0  (x,_nf);
      double P1   = PLLp1  (x,_nf);
      double c20 = xTCggH[_RCvar][k]->GetCoeff20();
      double c10 = xTCggH[_RCvar][k]->GetCoeff10();
      double c11 = xTCggH[_RCvar][k]->GetCoeff11();
      double b0 = beta0(_nf);
      //if(aux) res += as*as * (  c20*(P0sq-b0*P0) +   c10*P1);
      //else    res -= as*as * (2*c20*(P0sq-b0*P0) + 2*c10*P1 + c11*P0sq);
      if(aux) res -= as*as * (  c20*(P0sq-b0*P0) +   c10*P1);
      else    res += as*as * (2*c20*(P0sq-b0*P0) + 2*c10*P1 + c11*P0sq);
    }
    //if(matched_to_fixed_order == N3LO) {
    if(matched_to_fixed_order < N3LO) {
      double P0cb = PLLp0cb(x,_nf);
      double P0sq = PLLp0sq(x,_nf);
      double P0   = PLLp0  (x,_nf);
      double P1   = PLLp1  (x,_nf);
      //double P2   = PLLp2  (x,_nf);//,_RCvar); // here we don't have to vary it - we always use the same gammaLLp
      double P2   = PNLL2  (x,_nf);//,_RCvar); // here we don't have to vary it - we always use the same gammaLLp
      double P0P1 = PLLp0PLLp1(x,_nf);
      double c30 = xTCggH[_RCvar][k]->GetCoeff30();
      double c20 = xTCggH[_RCvar][k]->GetCoeff20();
      double c10 = xTCggH[_RCvar][k]->GetCoeff10();
      double c11 = xTCggH[_RCvar][k]->GetCoeff11();
      double c21 = xTCggH[_RCvar][k]->GetCoeff21();
      double b0 = beta0(_nf);
      double T = (_RCvar==1 ? 1 : 2);
      //if(aux) res += as*as*as * (  c30*(P0cb-3*b0*P0sq+2*b0*b0*P0) +   c20*(2*P0P1-T*b0*P1) +   c10*P2);
      //else    res -= as*as*as * (2*c30*(P0cb-3*b0*P0sq+2*b0*b0*P0) + 2*c20*(2*P0P1-T*b0*P1) + 2*c10*P2 + 2*c11*P0P1 + 2*c21*(P0cb-b0*P0sq));
      if(aux) res -= as*as*as * (  c30*(P0cb-3*b0*P0sq+2*b0*b0*P0) +   c20*(2*P0P1-T*b0*P1) +   c10*P2);
      else    res += as*as*as * (2*c30*(P0cb-3*b0*P0sq+2*b0*b0*P0) + 2*c20*(2*P0P1-T*b0*P1) + 2*c10*P2 + 2*c11*P0P1 + 2*c21*(P0cb-b0*P0sq));
    }
    return res;
  }


  double HELLxnf::deltaC2g  (double as, double x, Order matched_to_fixed_order) {
    return DeltaC(as, x, matched_to_fixed_order, "F2");
  }
  double HELLxnf::deltaC2q  (double as, double x, Order matched_to_fixed_order) {
    return deltaC2g(as, x, matched_to_fixed_order) * CF/CA;
  }

  double HELLxnf::deltaCLg  (double as, double x, Order matched_to_fixed_order) {
    return DeltaC(as, x, matched_to_fixed_order, "FL");
  }
  double HELLxnf::deltaCLq  (double as, double x, Order matched_to_fixed_order) {
    return deltaCLg(as, x, matched_to_fixed_order) * CF/CA;
  }

  double HELLxnf::deltaKhg (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return DeltaCm(as, x, matched_to_fixed_order, "Khg", m_Q_ratio);
  }
  double HELLxnf::deltaKhq (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaKhg(as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }

  double HELLxnf::deltaMC2g (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    if(4*m_Q_ratio*m_Q_ratio*x/(1-x)>1) return 0;
    if(m_Q_ratio < 0.002) return deltaC2g(as, x, matched_to_fixed_order)/_nf + deltaKhg(as, x, m_Q_ratio, matched_to_fixed_order);
    return DeltaCm(as, x, matched_to_fixed_order, "F2m", m_Q_ratio);
  }
  double HELLxnf::deltaMC2q (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaMC2g(as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }

  double HELLxnf::deltaMCLg (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    if(4*m_Q_ratio*m_Q_ratio*x/(1-x)>1) return 0;
    if(m_Q_ratio < 0.002) return deltaCLg(as, x, matched_to_fixed_order)/_nf;
    return DeltaCm(as, x, matched_to_fixed_order, "FLm", m_Q_ratio);
  }
  double HELLxnf::deltaMCLq (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaMCLg(as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }

  double HELLxnf::deltaMC2axg (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    if(4*m_Q_ratio*m_Q_ratio*x/(1-x)>1) return 0;
    if(m_Q_ratio < 0.002) return 0;
    return DeltaCm(as, x, matched_to_fixed_order, "F2axm", m_Q_ratio);
  }
  double HELLxnf::deltaMC2axq (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaMC2axg(as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }

  double HELLxnf::deltaMCLaxg (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    if(4*m_Q_ratio*m_Q_ratio*x/(1-x)>1) return 0;
    if(m_Q_ratio < 0.002) return 0;
    return DeltaCm(as, x, matched_to_fixed_order, "FLaxm", m_Q_ratio);
  }
  double HELLxnf::deltaMCLaxq (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaMCLaxg(as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }

  double HELLxnf::deltaMC2CCg (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    if(m_Q_ratio*m_Q_ratio*x/(1-x)>1) return 0;
    if(m_Q_ratio < 0.002) return deltaC2g(as, x, matched_to_fixed_order)/_nf + deltaKhg(as, x, m_Q_ratio, matched_to_fixed_order);
    return DeltaCm(as, x, matched_to_fixed_order, "F2CCm", m_Q_ratio);
  }
  double HELLxnf::deltaMC2CCq (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaMC2CCg(as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }

  double HELLxnf::deltaMCLCCg (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    if(m_Q_ratio*m_Q_ratio*x/(1-x)>1) return 0;
    if(m_Q_ratio < 0.002) return deltaCLg(as, x, matched_to_fixed_order)/_nf;
    return DeltaCm(as, x, matched_to_fixed_order, "FLCCm", m_Q_ratio);
  }
  double HELLxnf::deltaMCLCCq (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaMCLCCg(as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }

  double HELLxnf::deltaMC3CCg (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    if(m_Q_ratio*m_Q_ratio*x/(1-x)>1) return 0;
    if(m_Q_ratio < 0.002) return 0;
    return DeltaCm(as, x, matched_to_fixed_order, "F3CCm", m_Q_ratio);
  }
  double HELLxnf::deltaMC3CCq (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaMC3CCg(as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }


  double HELLxnf::deltaCggH   (double as, double x, double muFrat, Order matched_to_fixed_order) {
    return DeltaCggH(as, x, muFrat, matched_to_fixed_order, false);
  }
  double HELLxnf::deltaCggHaux(double as, double x, double muFrat, Order matched_to_fixed_order) {
    return DeltaCggH(as, x, muFrat, matched_to_fixed_order, true);
  }











  // **************************
  // ***  variable nf part  ***
  // **************************

  HELLx::HELLx(LogOrder order, string prepath) {
    for(int inf=0; inf<4; inf++)
      sxD[inf] = new HELLxnf(3+inf, order, prepath);
  }
  HELLx::~HELLx() {
    for(int inf=0; inf<4; inf++)
      delete sxD[inf];
  }

  void check_nf(int nf) {
    if(nf<3||nf>6) {
      cout << "\033[0;31m" << "HELLx: Non valid value of nf = " << nf << ". Allowed range nf=[3,6]." << "\033[0m" << endl;
      exit(234);
    }
  }

  HELLxnf* HELLx::GetHELLxnf(int nf) {
    check_nf(nf);
    return sxD[nf-3];
  }

    

  // Delta P Matrix
  sqmatrix<double> HELLx::DeltaP(int nf, double as, double x,  Order matched_to_fixed_order) {
    check_nf(nf);
    return sxD[nf-3]->DeltaP(as, x, matched_to_fixed_order);
  }


  double HELLx::deltaC2g  (int nf, double as, double x, Order matched_to_fixed_order) {
    check_nf(nf);
    return sxD[nf-3]->deltaC2g(as, x, matched_to_fixed_order);
  }
  double HELLx::deltaC2q  (int nf, double as, double x, Order matched_to_fixed_order) {
    return deltaC2g(nf, as, x, matched_to_fixed_order) * CF/CA;
  }

  double HELLx::deltaCLg  (int nf, double as, double x, Order matched_to_fixed_order) {
    check_nf(nf);
    return sxD[nf-3]->deltaCLg(as, x, matched_to_fixed_order);
  }
  double HELLx::deltaCLq  (int nf, double as, double x, Order matched_to_fixed_order) {
    return deltaCLg(nf, as, x, matched_to_fixed_order) * CF/CA;
  }

  double HELLx::deltaKhg (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    check_nf(nf);
    if(nf==6) { // Please check this!!!!
      cout << "HELLx: You requested matching function in the nf=6 scheme. Isn't it too much? Returning zero instead..." << endl;
      return 0;
    }
    return sxD[nf-3]->deltaKhg(as, x, m_Q_ratio, matched_to_fixed_order);
  }
  double HELLx::deltaKhq (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaKhg(nf, as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }

  double HELLx::deltaMC2g (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    check_nf(nf);
    if(nf==6) {
      cout << "HELLx: You requested massive CFs in the nf=6 scheme. Isn't it too much? Returning massless coefficient functions instead..." << endl;
      return deltaC2g(nf, as, x, matched_to_fixed_order);
    }
    return sxD[nf-3]->deltaMC2g(as, x, m_Q_ratio, matched_to_fixed_order);
  }
  double HELLx::deltaMC2q (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaMC2g(nf, as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }

  double HELLx::deltaMCLg (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    check_nf(nf);
    if(nf==6) {
      cout << "HELLx: You requested massive CFs in the nf=6 scheme. Isn't it too much? Returning massless coefficient functions instead..." << endl;
      return deltaCLg(nf, as, x, matched_to_fixed_order);
    }
    return sxD[nf-3]->deltaMCLg(as, x, m_Q_ratio, matched_to_fixed_order);
  }
  double HELLx::deltaMCLq (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaMCLg(nf, as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }

  double HELLx::deltaMC2axg (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    check_nf(nf);
    if(nf==6) {
      cout << "HELLx: You requested massive CFs in the nf=6 scheme. Isn't it too much? Returning zero instead..." << endl;
      return 0;
    }
    return sxD[nf-3]->deltaMC2axg(as, x, m_Q_ratio, matched_to_fixed_order);
  }
  double HELLx::deltaMC2axq (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaMC2axg(nf, as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }

  double HELLx::deltaMCLaxg (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    check_nf(nf);
    if(nf==6) {
      cout << "HELLx: You requested massive CFs in the nf=6 scheme. Isn't it too much? Returning zero instead..." << endl;
      return 0;
    }
    return sxD[nf-3]->deltaMCLaxg(as, x, m_Q_ratio, matched_to_fixed_order);
  }
  double HELLx::deltaMCLaxq (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaMCLaxg(nf, as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }

  double HELLx::deltaMC2CCg (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    check_nf(nf);
    if(nf==6) {
      cout << "HELLx: You requested massive CFs in the nf=6 scheme. Isn't it too much? Returning massless coefficient functions instead..." << endl;
      return deltaC2g(nf, as, x, matched_to_fixed_order);
    }
    return sxD[nf-3]->deltaMC2CCg(as, x, m_Q_ratio, matched_to_fixed_order);
  }
  double HELLx::deltaMC2CCq (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaMC2CCg(nf, as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }

  double HELLx::deltaMCLCCg (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    check_nf(nf);
    if(nf==6) {
      cout << "HELLx: You requested massive CFs in the nf=6 scheme. Isn't it too much? Returning massless coefficient functions instead..." << endl;
      return deltaCLg(nf, as, x, matched_to_fixed_order);
    }
    return sxD[nf-3]->deltaMCLCCg(as, x, m_Q_ratio, matched_to_fixed_order);
  }
  double HELLx::deltaMCLCCq (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaMCLCCg(nf, as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }

    double HELLx::deltaMC3CCg (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    check_nf(nf);
    if(nf==6) {
      cout << "HELLx: You requested massive CFs in the nf=6 scheme. Isn't it too much? Returning zero instead..." << endl;
      return 0;
    }
    return sxD[nf-3]->deltaMC3CCg(as, x, m_Q_ratio, matched_to_fixed_order);
  }
  double HELLx::deltaMC3CCq (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaMC3CCg(nf, as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }












};

