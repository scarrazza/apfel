#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>

#include "HELL/include/hell-x.hh"
#include "HELL/include/expansionSFs.hh"
#include "HELL/include/math/special_functions.hh"


using namespace std;



namespace HELLx {

  bool _fast_ = true;
  bool a10woRC = false;
  bool fullySymm = false;
  bool shiftMmin = true;
  
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
  //inline double c1NLO(int nf) { return -15.04043249 - 0.2076545137*nf; }
  //inline double k1NLO(int nf) { return -574.496771  - 0.3724353405*nf; }
  //inline double m1NLO(int nf) { return 0.3099881585 - 0.01878716112*nf; }
  inline double c1NLO(int nf) { return -11.696833425 - 0.4102968810*nf; }
  inline double k1NLO(int nf) { return -494.250393369 - 5.23585215538*nf; }
  inline double m1NLO(int nf) { return -0.1492429211 + 0.00904502552*nf; }
  //
  inline double a11()       { return CA/M_PI; }
  inline double a10(int nf) { return -(11.*CA + 2.*nf*(1.-2.*CF/CA))/12./M_PI + (a10woRC ? beta0(nf) : 0 ); }
  inline double a21(int nf) { return nf*(26*CF-23*CA)/36/M_PI/M_PI; }
  //
  unsigned int factorials[] = {1,1,2,6,24,120,720};
  unsigned int factorial(unsigned int k) {
    if(k<7) return factorials[k];
    double res = 1;
    for(int j=1; j<=(int) k; j++) res *= j;
    return res;
  }
  double binomial(unsigned int k, unsigned int j) {
    return factorial(k)/(1.*factorial(k-j)*factorial(j));
  }
  double Pole(double x, int k, int j, int p=0) { // inverse Mellin of 1/N^k/(1+N)^j/(2+N)^p
    if(k==0 && j==0 && p>0) return pow(-log(x),p-1)/factorial(p-1)*x;
    if(p==0 && k==0 && j>0) return pow(-log(x),j-1)/factorial(j-1);
    if(p==0 && j==0 && k>0) return pow(-log(x),k-1)/factorial(k-1)/x;
    if(p==0 && k==1 && j==1) return 1/x-1;
    if(p==0 && k==2 && j==1) return Pole(x,2,0) - Pole(x,1,0) + Pole(x,0,1);
    if(p==0 && k==1 && j==2) return Pole(x,1,0) - Pole(x,0,2) - Pole(x,0,1);
    cout << "HELLx warning: this inverse mellin is not implemented: 1/N^"<<k<<"/(1+N)^"<<j<<"/(2+N)^"<<p << endl;
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




  // expansion of the LLp/NLL
  double Paux0(double x, int nf) {
    double a1 = a11(), a0=a10(nf) - (a10woRC ? beta0(nf) : 0 );
    double res = a1/x + a0;
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k) * pow(1-sqrt(x),j);
  }
  double Paux0sq(double x, int nf) {
    double a1 = a11(), a0=a10(nf) - (a10woRC ? beta0(nf) : 0 );
    double res = a1*a1*Pole(x,2,0) + a0*a0*Pole(x,0,2) + 2*a1*a0*Pole(x,1,1);
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k) * pow(1-sqrt(x),j);
  }
  double Paux0cb(double x, int nf) {
    double a1 = a11(), a0=a10(nf) - (a10woRC ? beta0(nf) : 0 );
    double res = a1*a1*a1*Pole(x,3,0) + a0*a0*a0*Pole(x,0,3) + 3*a1*a1*a0*Pole(x,2,1) + 3*a1*a0*a0*Pole(x,1,2);
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k) * pow(1-sqrt(x),j);
  }
  double Paux1(double x, int nf, bool useLLp) {
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
  double Paux0Paux1(double x, int nf, bool useLLp) {
    double a1 = a11(), a0=a10(nf) - (a10woRC ? beta0(nf) : 0 );
    double res = 0;
    if(useLLp) {
      res = beta0(nf)*(3*k0/32.-c0) * ( a1*(Pole(x,2,0)-4*Pole(x,0,2)) + a0*(Pole(x,1,0)+4*Pole(x,0,3)-4*Pole(x,0,2)-Pole(x,0,1)) );
    } else {
      res = a21(nf) * ( a1*Pole(x,2,0) + (a0-2*a1)/x -2*a0*Pole(x,0,2) + 2*a1-a0 );
    }
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k) * pow(1-sqrt(x),j);
  }
  double Paux2(double x, int nf, bool useLLp, int RCvar) {
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
      int T = (RCvar==1 ? 1 : 2);
      double  c1 = c1NLO(nf);
      double  k1 = k1NLO(nf);
      double  m1 = m1NLO(nf);
      double dc1 = c1-c1LO(nf);
      double dk1 = k1-k1LO(nf);
      double dm1 = m1;
      if(shiftMmin) m1 += c0/2.;
      double lam2 = b0*(48*c0*k0-3*k0*k0-256*c0*c0)/256. + (16*c0*(dk1+6*k0*dm1)+k0*(16*dc1-3*dk1-15*k0*dm1))/512.;
      double lam1 = b0*(3*T*k1/32 -T*c1 + (b0+6*m1)*k0/16.);
      res += lam1*Pole(x,1,0) + lam2*Pole(x,2,0) - 4*(1+log(x))*(lam1+lam2);
    } else {
      double chi01 = a1;
      double chi02 = a0*a1;
      double chi11 = a21(nf);
      double A = chi01*chi01*(5*ZETA3/2.-1./24.) + chi02*(53./18.-ZETA2) - chi11;
      //double A = ( CA*CA*(-74./27. +11./12.*ZETA2 +5./2.*ZETA3) +nf*((4*CA+7*CF)/27. +(CA-2*CF)*ZETA2/6.) )/M_PI/M_PI;
      res = (A*a1 + a1*chi02 + a1*chi11) * Pole(x,2,0)
	+ (-4*A*a1 - 3*a1*chi02 - 3*a1*chi11) * Pole(x,1,0)
	+ 4*(a0 + a1)*(2*A + chi02 + chi11) * Pole(x,0,4)
	- 2*(a0 + a1)*(4*A + chi02 + chi11) * Pole(x,0,3)
	+ 2*(A*a0 + 3*A*a1 + a1*chi02 + a1*chi11) * Pole(x,0,2)
	+ (4*A*a1 + 3*a1*chi02 + 3*a1*chi11) * Pole(x,0,1)
	;
      res += chi01 * ( a1*a1*  exactPoly(x,3,0) 
		       +a1*a0* exactPoly(x,2,0)
		       +4*pow(a0 + a1,2)*              exactPoly(x,0,3)
		       -2*(a0*a0 + 4*a0*a1 + 3*a1*a1)* exactPoly(x,0,2)
		       );
      if(!a10woRC && !fullySymm) res += b0 * chi01 * (a1*exactPoly(x,2,0) + 2*(a1+a0)*exactPoly(x,0,2) + a1*4*(Pole(x,1,0)-Pole(x,0,1)-Pole(x,0,2)) +2*(a1+a0)*(Pole(x,0,3)-Pole(x,0,4)) );
      // add RC
      res += b0*b0*k0/16* (1/x-4*(1+log(x)));
    }
    // result
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k)*pow(1-sqrt(x),j);
  }






  // expansion of the LL
  double PLL1(double x, int nf) { return Paux1(x,nf,true); }
  double PLL2(double x, int nf, int RCvar) {
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
    int T = (RCvar==1 ? 1 : 2);
    double b0 = beta0(nf);
    double c1 = c1LO(nf);
    double k1 = k1LO(nf);
    double m1 = 0.5;
    if(shiftMmin) m1 += c0/2.;
    double lam2 = b0*(48*c0*k0-3*k0*k0-256*c0*c0)/256.;
    double lam1 = b0*(3*T*k1/32 -T*c1 + (b0+6*m1)*k0/16.);
    res += lam1*Pole(x,1,0) + lam2*Pole(x,2,0) - 4*(1+log(x))*(lam1+lam2);
    // result
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k)*pow(1-sqrt(x),j);
  }






  // expansion of the NLL
  double PNLL2(double x, int nf) {  // this is independent of var
    double res = 0;
    double b0 = beta0(nf);
    double a1 = a11(), a0=a10(nf);
    double chi01 = a1;
    double chi02 = a0*a1;
    double chi11 = a21(nf);
    double A = chi01*chi01*(5*ZETA3/2.-1./24.) + chi02*(53./18.-ZETA2) - chi11;
    //double A = ( CA*CA*(-74./27. +11./12.*ZETA2 +5./2.*ZETA3) +nf*((4*CA+7*CF)/27. +(CA-2*CF)*ZETA2/6.) )/M_PI/M_PI;
    res = (A*a1 + a1*chi02 + a1*chi11) * Pole(x,2,0)
      + (-4*A*a1 - 3*a1*chi02 - 3*a1*chi11) * Pole(x,1,0)
      + 4*(a0 + a1)*(2*A + chi02 + chi11) * Pole(x,0,4)
      - 2*(a0 + a1)*(4*A + chi02 + chi11) * Pole(x,0,3)
      + 2*(A*a0 + 3*A*a1 + a1*chi02 + a1*chi11) * Pole(x,0,2)
      + (4*A*a1 + 3*a1*chi02 + 3*a1*chi11) * Pole(x,0,1)
      ;
    if(_fast_) {
      res += CA/M_PI* ( a1*a1*  Poly(x,3,0) 
			+a1*a0* Poly(x,2,0)
			+4*pow(a0 + a1,2)*              Poly(x,0,3)
			-2*(a0*a0 + 4*a0*a1 + 3*a1*a1)* Poly(x,0,2)
			);
      if(!a10woRC && !fullySymm) res += b0 * chi01 * (a1*Poly(x,2,0) + 2*(a1+a0)*Poly(x,0,2) + a1*4*(Pole(x,1,0)-Pole(x,0,1)-Pole(x,0,2)) +2*(a1+a0)*(Pole(x,0,3)-Pole(x,0,4)) );
    } else {
      res += CA/M_PI* ( a1*a1*  exactPoly(x,3,0) 
			+a1*a0* exactPoly(x,2,0)
			+4*pow(a0 + a1,2)*              exactPoly(x,0,3)
			-2*(a0*a0 + 4*a0*a1 + 3*a1*a1)* exactPoly(x,0,2)
			);
      if(!a10woRC && !fullySymm) res += b0 * chi01 * (a1*exactPoly(x,2,0) + 2*(a1+a0)*exactPoly(x,0,2) + a1*4*(Pole(x,1,0)-Pole(x,0,1)-Pole(x,0,2)) +2*(a1+a0)*(Pole(x,0,3)-Pole(x,0,4)) );
    }
    // add RC
    res += b0*b0*k0/16* (1/x-4*(1+log(x)));
    // result
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k)*pow(1-sqrt(x),j);
  }
  // coefficients of the pole expansion in N=0,-1,-2
  double P3exp000[4] = {-29.15790, 24.5597755, -15.90479605, 1.999127764};
  double P3exp001[4] = {-23.669696, 6.77716360, 0, 0};
  double P3exp010[4] = {-1.8220376, 0.099304430, 0.504630968, 0};
  double P3exp011[4] = {2.6183086, -0.1730105325, 0, 0};
  double P3exp020[4] = {-0.032450437, 0.0144922214, 0, 0};
  double P3exp021[4] = {-0.056547647, -0.01440767607, 0, 0};
  double P3exp030[4] = {0.00100208876, 0, 0, 0};
  double P3exp031[4] = {-0.00092101170, 0, 0, 0};
  //
  double P3exp100[7] = {33.8085, 29.2302, -6.69277937, 1.834031285, 1.915741435, -0.1174171720, 0.00384974335};
  double P3exp101[7] = {67.5701, -67.570, 0, 0, 0, 0, 0};
  double P3exp110[7] = {1.22928,1.23230,-1.017298598,1.235048578,-0.1769830668,0.01611188882,-0.000855498521};
  double P3exp111[7] = {-9.7812, 9.7812, 0, 0, 0, 0, 0};
  double P3exp120[7] = {0.104026,-0.151810,0.0969643030,-0.0736316691,0.00319322032,-0.000454153536,0.0000633702608};
  double P3exp121[7] = {0.283821, -0.283821, 0, 0, 0, 0, 0};
  double P3exp130[7] = {-0.00405723,0.0042226,-0.000542130102,0.0003528635814,-0.0000440723213,-7.04114009e-6,-1.564697798e-6};
  double P3exp131[7] = {0.00368405, -0.0036840, 0, 0, 0, 0, 0};
  //
  double P3exp200[3] = {-11.10823839, -4.510900857, 0.1212669154};
  double P3exp210[3] = {0.776501741, 0.3370529857, 0.01458624979};
  double P3exp220[3] = {-0.01026240061, -0.002119242617, -0.000744600565};
  double P3exp230[3] = {0.0000542637196, 0.00001080423830, 6.10232141e-6};
  //
  double PNLL3(double x, int nf, int RCvar) {
    double res = 0;
    int T = (RCvar==1 ? 1 : 2);
    for(int k=0; k<4; k++) {
      res += P3exp000[k] * Pole(x,k+1,0,0);
      res += P3exp001[k] * Pole(x,k+1,0,0) * T;
      res += P3exp010[k] * Pole(x,k+1,0,0) * nf;
      res += P3exp011[k] * Pole(x,k+1,0,0) * nf * T;
      res += P3exp020[k] * Pole(x,k+1,0,0) * nf*nf;
      res += P3exp021[k] * Pole(x,k+1,0,0) * nf*nf * T;
      res += P3exp030[k] * Pole(x,k+1,0,0) * nf*nf*nf;
      res += P3exp031[k] * Pole(x,k+1,0,0) * nf*nf*nf * T;
    }
    for(int k=0; k<7; k++) {
      res += P3exp100[k] * Pole(x,0,k+1,0);
      res += P3exp101[k] * Pole(x,0,k+1,0) * T;
      res += P3exp110[k] * Pole(x,0,k+1,0) * nf;
      res += P3exp111[k] * Pole(x,0,k+1,0) * nf * T;
      res += P3exp120[k] * Pole(x,0,k+1,0) * nf*nf;
      res += P3exp121[k] * Pole(x,0,k+1,0) * nf*nf * T;
      res += P3exp130[k] * Pole(x,0,k+1,0) * nf*nf*nf;
      res += P3exp131[k] * Pole(x,0,k+1,0) * nf*nf*nf * T;
    }
    for(int k=0; k<3; k++) {
      res += P3exp100[k] * Pole(x,0,0,k+1);
      res += P3exp110[k] * Pole(x,0,0,k+1) * nf;
      res += P3exp120[k] * Pole(x,0,0,k+1) * nf*nf;
      res += P3exp130[k] * Pole(x,0,0,k+1) * nf*nf*nf;
    }
    // result
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k)*pow(1-sqrt(x),j);
  }







  // anomalous dimensions, needed for momcons
  // LLp/NLL
  double gamma0aux(double N, int nf) {
    double a1 = a11(), a0=a10(nf);
    return a1/N + a0/(N+1);
  }
  double gamma1aux(double N, int nf, bool useLLp) {
    double b0 = beta0(nf);
    return (useLLp ? b0*(3*k0/32.-c0) * (1/N - fmom(N)) : a21(nf) * (1/N-2/(N+1)) );
  }
  double gamma2aux(double N, int nf, bool useLLp, bool fast, int RCvar) {
    double b0 = beta0(nf);
    double res = 0;
    if(useLLp) {
      double a1 = a11(), a0=a10(nf);
      double chi01 = a1;
      double chi02 = a1*a0;
      double gamma0  = a1/N + a0 - 2*(a1+a0)*N/(N+1);
      double gamma0p = -a1/N/N - 2*(a1+a0)/(N+1)/(N+1);
      res  = -gamma0p * ( chi02/pow(N+1.,2) - chi02/4.*fmom(N) + chi01*gamma0*(dpsi(N+1)-ZETA2) );
      // add RC
      int T = (RCvar==1 ? 1 : 2);
      double  c1 = c1NLO(nf);
      double  k1 = k1NLO(nf);
      double  m1 = m1NLO(nf);
      double dc1 = c1-c1LO(nf);
      double dk1 = k1-k1LO(nf);
      double dm1 = m1;
      if(shiftMmin) m1 += c0/2.;
      double lam2 = b0*(48*c0*k0-3*k0*k0-256*c0*c0)/256. + (16*c0*(dk1+6*k0*dm1)+k0*(16*dc1-3*dk1-15*k0*dm1))/512.;
      double lam1 = b0*(3*T*k1/32 -T*c1 + (b0+6*m1)*k0/16.);
      res += lam2/N/N + lam1/N - (lam1+lam2)*fmom(N);
    }
    else { // NLL
      double a1 = a11(), a0=a10(nf);
      double chi01 = a1;
      double chi02 = a1*a0;
      double chi11 = a21(nf);
      double A = chi01*chi01*(5*ZETA3/2.-1./24.) + chi02*(53./18.-ZETA2) - chi11;
      //double A = ( CA*CA*(-74./27. +11./12.*ZETA2 +5./2.*ZETA3) +nf*((4*CA+7*CF)/27. +(CA-2*CF)*ZETA2/6.) )/M_PI/M_PI;
      double gamma0p = -a1/N/N - 2*(a1+a0)/(N+1)/(N+1);
      res  = -gamma0p * ( A + chi11/(N+1.) + chi02/pow(N+1.,2) - (A+chi11/2.+chi02/4.)*fmom(N) );
      if(fast) {
	res += chi01* ( a1*a1*  ( 3*ZETA4/N -2*ZETA3/N/N -3/(N+1) -1/(N+1)/(N+1) )
			+a1*a0* (-2*ZETA3/N +2/(N+1) +1/(N+1)/(N+1) )
			+4*pow(a0 + a1,2)*              ( 3*ZETA4/(N+1) -2*ZETA3/(N+1)/(N+1) +1/pow(N+1,5) )
			-2*(a0*a0 + 4*a0*a1 + 3*a1*a1)* ( -2*ZETA3/(N+1) +1/pow(N+1,4) )
			);
	if(!a10woRC && !fullySymm) res -= b0*chi01* ( -a1*(-2*ZETA3/N+2/(N+1)+1/(N+1)/(N+1)) - 2*(a1+a0)*(-2*ZETA3/(N+1)+1/pow(N+1,4)) + gamma0p*fmom(N) );
      } else {
	double gamma0  = a1/N + a0 - 2*(a1+a0)*N/(N+1);
	res -= chi01 * gamma0 * gamma0p * (dpsi(N+1)-ZETA2);
	if(!a10woRC && !fullySymm) res -= gamma0p * b0 * chi01 * (dpsi(N+1)-ZETA2 + fmom(N));
      }
      // add RC
      res += b0*b0*k0/16.* (1/N - fmom(N));
    }
    return res;
  }
  // gg LL
  double gammagg1LL(double N, int nf) {
    return gamma1aux(N,nf,true); //beta0(nf)* (3*k0/32.-c0) * (1./N - fmom(N));
  }
  double gammagg2LL(double N, int nf, int RCvar) {
    double res = 0;
    double a1 = a11(), a0=a10(nf);
    double chi02 = a1*a0;
    double gamma0  = a1/N + a0 - 2*(a1+a0)*N/(N+1);
    double gamma0p = -a1/N/N - 2*(a1+a0)/(N+1)/(N+1);
    res  = -gamma0p * ( chi02/pow(N+1.,2) - chi02/4*fmom(N) + a1*gamma0*(dpsi(N+1)-ZETA2) );
    // add RC
    int T = (RCvar==1 ? 1 : 2);
    double b0 = beta0(nf);
    double rho2 = b0*(48*c0*k0-3*k0*k0-256*c0*c0)/256.;
    double rho1 = b0*(3*T*k1LO(nf)/32 -T*c1LO(nf) + b0*k0/16);
    res += rho2/N/N + rho1/N - (rho2+rho1)*fmom(N);
    // result
    return res;
  }
  // gg NLL
  double gammagg2NLL(double N, int nf, bool useLLp) { // this is independent of var
    double res = gamma2aux(N,nf,false,_fast_,0);
    // add Pqg
    double b0 = beta0(nf);
    double gamma0LL = gamma0aux(N,nf); //a1/N + a0/(N+1);
    double gamma1LL = gamma1aux(N,nf,useLLp); //(useLLp ? b0*(3*k0/32.-c0) * (1/N - fmom(N)) : a21(nf) * (1/N-2/(N+1)) );
    double h1 = 5./3., h2 = 14./9.;
    res += (1-CF/CA) * (h2*gamma0LL*(gamma0LL-b0) + h1*gamma1LL) * nf/3./M_PI;
    // result
    return res;
  }
  double gammagg3NLL(double N, int nf, bool useLLp, int RCvar) {
    double res = 0;
    int T = (RCvar==1 ? 1 : 2);
    for(int k=0; k<4; k++) {
      res += P3exp000[k] * pow(N,-k-1);
      res += P3exp001[k] * pow(N,-k-1) * T;
      res += P3exp010[k] * pow(N,-k-1) * nf;
      res += P3exp011[k] * pow(N,-k-1) * nf * T;
      res += P3exp020[k] * pow(N,-k-1) * nf*nf;
      res += P3exp021[k] * pow(N,-k-1) * nf*nf * T;
      res += P3exp030[k] * pow(N,-k-1) * nf*nf*nf;
      res += P3exp031[k] * pow(N,-k-1) * nf*nf*nf * T;
    }
    for(int k=0; k<7; k++) {
      res += P3exp100[k] * pow(N+1,-k-1);
      res += P3exp101[k] * pow(N+1,-k-1) * T;
      res += P3exp110[k] * pow(N+1,-k-1) * nf;
      res += P3exp111[k] * pow(N+1,-k-1) * nf * T;
      res += P3exp120[k] * pow(N+1,-k-1) * nf*nf;
      res += P3exp121[k] * pow(N+1,-k-1) * nf*nf * T;
      res += P3exp130[k] * pow(N+1,-k-1) * nf*nf*nf;
      res += P3exp131[k] * pow(N+1,-k-1) * nf*nf*nf * T;
    }
    for(int k=0; k<3; k++) {
      res += P3exp100[k] * pow(N+2,-k-1);
      res += P3exp110[k] * pow(N+2,-k-1) * nf;
      res += P3exp120[k] * pow(N+2,-k-1) * nf*nf;
      res += P3exp130[k] * pow(N+2,-k-1) * nf*nf*nf;
    }
    // add Pqg
    int T2 = (RCvar==2 ? 1 : 2);  // Pqg variation
    double b0 = beta0(nf);
    double gamma0LL = gamma0aux(N,nf);
    double gamma1LL = gamma1aux(N,nf,useLLp);
    double gamma2LL = gamma2aux(N,nf,useLLp,_fast_,RCvar);  // here we don;t have to vary it, as in Pqg we always use the same gammaLLp/NLL
    double h1 = 5./3., h2 = 14./9., h3 = 2*(41./81.+ZETA3);
    res += (1-CF/CA) * (h3*gamma0LL*(gamma0LL-b0)*(gamma0LL-2*b0) + h2*gamma1LL*(2*gamma0LL-T2*b0) + h1*gamma2LL) * nf/3./M_PI;
    // result
    return res;
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
  double mcPgg2LL(double x, int nf, int RCvar) {
    double momcons = 0;
    int k=_damping, j=_dampingsqrt;
    if(_damping==2 && _dampingsqrt==4) {
      double coeffs[9] = {1, -4, 4, 4, -10, 4, 4, -4, 1};
      for(int n=0; n<=8; n++) {
	momcons += coeffs[n] * gammagg2LL(1.+n/2.,nf,RCvar);
      }
    } else {
      for(int n=0; n<=k; n++) {
	for(int m=0; m<=j; m++) {
	  momcons += binomial(k,n) * binomial(j,m) * pow(-1.,n+m) * gammagg2LL(1.+n+m/2.,nf,RCvar);
	}
      }
    }
    momcons /= d1(k,j);
    return momcons * pow(1-x,k)*pow(1-sqrt(x),j);
  }
  double mcPgg2NLL(double x, int nf, bool useLLp) {
    double momcons = 0;
    int k=_damping, j=_dampingsqrt;
    if(_damping==2 && _dampingsqrt==4) {
      double coeffs[9] = {1, -4, 4, 4, -10, 4, 4, -4, 1};
      for(int n=0; n<=8; n++) {
	momcons += coeffs[n] * gammagg2NLL(1.+n/2.,nf,useLLp);
      }
    } else {
      for(int n=0; n<=k; n++) {
	for(int m=0; m<=j; m++) {
	  momcons += binomial(k,n) * binomial(j,m) * pow(-1.,n+m) * gammagg2NLL(1.+n+m/2.,nf,useLLp);
	}
      }
    }
    momcons /= d1(k,j);
    return momcons * pow(1-x,k)*pow(1-sqrt(x),j);
  }
  double mcPgg3NLL(double x, int nf, bool useLLp, int RCvar) {
    double momcons = 0;
    int k=_damping, j=_dampingsqrt;
    if(_damping==2 && _dampingsqrt==4) {
      double coeffs[9] = {1, -4, 4, 4, -10, 4, 4, -4, 1};
      for(int n=0; n<=8; n++) {
	momcons += coeffs[n] * gammagg3NLL(1.+n/2.,nf,useLLp,RCvar);
      }
    } else {
      for(int n=0; n<=k; n++) {
	for(int m=0; m<=j; m++) {
	  momcons += binomial(k,n) * binomial(j,m) * pow(-1.,n+m) * gammagg3NLL(1.+n+m/2.,nf,useLLp,RCvar);
	}
      }
    }
    momcons /= d1(k,j);
    return momcons * pow(1-x,k)*pow(1-sqrt(x),j);
  }









  // Pqg
  double Pqg2(double x, int nf, bool useLLp) {
    double h1 = 5./3., h2 = 14./9.;
    double b0 = beta0(nf);
    double P0sq = Paux0sq(x,nf);
    double P0   = Paux0  (x,nf);
    double P1   = Paux1  (x,nf,useLLp);
    double res = h2*(P0sq-b0*P0) + h1*P1;
    return res*nf/3./M_PI;
  }
  double Pqg3(double x, int nf, bool useLLp, int RCvar) {
    double h1 = 5./3., h2 = 14./9., h3 = 2*(41./81.+ZETA3);
    double b0 = beta0(nf);
    double P0cb = Paux0cb(x,nf);
    double P0sq = Paux0sq(x,nf);
    double P0   = Paux0  (x,nf);
    double P1   = Paux1  (x,nf,useLLp);
    double P2   = Paux2  (x,nf,useLLp,RCvar);  // here we don;t have to vary it, as in Pqg we always use the same gammaLLp/NLL
    double P0P1 = Paux0Paux1(x,nf,useLLp);
    double T = (RCvar==2 ? 1 : 2);  // Pqg variation
    double res = h3*(P0cb-3*b0*P0sq+2*b0*b0*P0) + h2*(2*P0P1-T*b0*P1) + h1*P2;
    return res*nf/3./M_PI;
  }



};

