#include "HELL/include/gammaNLO.hh"

const double Nc = 3.;
const double CA = Nc;
const double TF = 0.5;
const double CF = (Nc*Nc-1.)/2./Nc;

//const double Pi2 = M_PI*M_PI;
const double fourPi2 = 16.*M_PI*M_PI;
//const double fourPi3 = 64.*M_PI*M_PI*M_PI;

const double ZETA2  =  1.64493406684822643647;
const double ZETA3  =  1.2020569031595942855;
//const double ZETA4  =  1.082323233711138191516;
//const double ZETA5  =  1.036927755143369926331;
const double EulerGamma = 0.5772156649015328606065120900824;


double ABS(double a){ return (a >0 ? a: -a); }

dcomplex psi(dcomplex Z){
  dcomplex SUB = 0. ;
  dcomplex ZZ = Z;
  if(ABS(imag(ZZ))<10.) { // if too close to the real axis...
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

dcomplex dpsi(dcomplex Z, int M){
  int K1, K2;
  dcomplex SUB = 0. , SUBM;
  dcomplex ZZ = Z;
  if(ABS(imag(ZZ))<10.) { // if too close to the real axis...
  label1:
    SUBM = -1./ZZ;
    for(K1=1; K1<=M; K1++) {
      SUBM = - SUBM * (double)K1 / ZZ;
    }
    if(real(ZZ)<10.) { // ...use recurrence relation to push real(z) large enough
      SUB = SUB + SUBM;
      ZZ = ZZ + 1.;
      goto label1;
    }
  }
  // Expansion coefficients for the first derivative
  double A1 =  1.;
  double A2 =  1./2.;
  double A3 =  1./6.;
  double A4 = -1./30.;
  double A5 =  1./42.;
  double A6 = -1./30.;
  double A7 =  5./66.;
  
  // Expansion coefficients for the higher derivatives
  if(M>1) {
    for(K2=2;K2<=M;K2++){
      A1 = A1 * (1.*K2-1.);
      A2 = A2 *  1.*K2;
      A3 = A3 * (1.*K2+1.);
      A4 = A4 * (1.*K2+3.);
      A5 = A5 * (1.*K2+5.);
      A6 = A6 * (1.*K2+7.);
      A7 = A7 * (1.*K2+9.);
    }
  }
  dcomplex RZ = 1./ ZZ;
  dcomplex DZ = RZ * RZ;
  // SUB + asympt expansion (Abramowitz, Stengun, 6.4.11)
  return SUB + pow(-1.,M+1) * pow(RZ,M) * ( A1 + RZ * (A2 + RZ * (A3 + DZ * (A4 + DZ * (A5 + DZ * (A6 + A7 * DZ ))))) );
}





dcomplex gamma0qq(dcomplex N) {
  return CF/2./M_PI * ( 1./N - 1./(N+1.) - 2.*psi(N+1.) - 2.*EulerGamma + 3./2. );
}
dcomplex gamma0qg(dcomplex N, double nf) {
  return nf/2./M_PI * ( 2.+N+N*N ) / N/(N+1.)/(N+2.);
}
dcomplex gamma0gq(dcomplex N) {
  return CF/2./M_PI * ( 2.+N+N*N ) / N/(N*N-1.);
}
dcomplex gamma0gg(dcomplex N, double nf) {
  return CA/M_PI * ( 11./12. + 1./(N-1.) - 1./N + 1./(N+1.) - 1./(N+2.) - psi(N+1.) - EulerGamma ) - nf/6./M_PI;
}






// Analytic continuations of the occuring sums
void sums(dcomplex N, gamma1sums &g1s) {
  //
  g1s.S1  = EulerGamma + psi(N+1.);
  g1s.S2  = ZETA2 - dpsi(N+1.,1);
  //    
  g1s.NS  = N * N;
  g1s.NT  = g1s.NS * N;
  g1s.NFO = g1s.NT * N;
  g1s.NFI = g1s.NFO * N;
  g1s.NSI = g1s.NFI * N;
  g1s.NSE = g1s.NSI * N;
  g1s.NE  = g1s.NSE * N;
  g1s.NN  = g1s.NE * N;
  // 
  g1s.NM  = N - 1.;
  g1s.N1  = N + 1.;
  g1s.N2  = N + 2.;
  g1s.NMS = g1s.NM * g1s.NM;
  g1s.N1S = g1s.N1 * g1s.N1;
  g1s.N1T = g1s.N1S * g1s.N1;
  g1s.N2S = g1s.N2 * g1s.N2;
  g1s.N2T = g1s.N2S * g1s.N2;
  //
  g1s.N3  = N + 3.;
  g1s.N4  = N + 4.;
  g1s.N5  = N + 5.;
  g1s.N6  = N + 6.;
  g1s.S11 = g1s.S1  + 1./g1s.N1;
  g1s.S12 = g1s.S11 + 1./g1s.N2;
  g1s.S13 = g1s.S12 + 1./g1s.N3;
  g1s.S14 = g1s.S13 + 1./g1s.N4;
  g1s.S15 = g1s.S14 + 1./g1s.N5;
  g1s.S16 = g1s.S15 + 1./g1s.N6;
  g1s.SPMOM = ( 1.0000 * (ZETA2 - g1s.S1 / N ) / N  -
		0.9992 * (ZETA2 - g1s.S11/ g1s.N1) / g1s.N1 +
		0.9851 * (ZETA2 - g1s.S12/ g1s.N2) / g1s.N2 -
		0.9005 * (ZETA2 - g1s.S13/ g1s.N3) / g1s.N3 +
		0.6621 * (ZETA2 - g1s.S14/ g1s.N4) / g1s.N4 -
		0.3174 * (ZETA2 - g1s.S15/ g1s.N5) / g1s.N5 +
		0.0699 * (ZETA2 - g1s.S16/ g1s.N6) / g1s.N6 );
  //
  g1s.SLC = - 5./8. * ZETA3;
  g1s.SLV = - ZETA2/2.* (psi(g1s.N1/2.) - psi(N/2.)) + g1s.S1/g1s.NS + g1s.SPMOM;
  g1s.SSCHLM = g1s.SLC - g1s.SLV;
  g1s.SSTR2M = ZETA2 - dpsi (g1s.N1/2.,1);
  g1s.SSTR3M = 0.5 * dpsi (g1s.N1/2.,2) + ZETA3;
  //
  g1s.SSCHLP = g1s.SLC + g1s.SLV;
  g1s.SSTR2P = ZETA2 - dpsi(g1s.N2/2.,1);
  g1s.SSTR3P = 0.5 * dpsi(g1s.N2/2.,2) + ZETA3;
  //
  return;
}


// P1NS
void P1NSplus(dcomplex N, gamma1sums g1s, dcomplex &PNPA, dcomplex &PNSB, dcomplex &PNSC) {
  PNPA = ( 16.* g1s.S1 * (2.* N + 1.) / (g1s.NS * g1s.N1S) +
	   16.* (2.* g1s.S1 - 1./(N * g1s.N1)) * ( g1s.S2 - g1s.SSTR2P ) +
	   64.* g1s.SSCHLP + 24.* g1s.S2 - 3. - 8.* g1s.SSTR3P -
	   8.* (3.* g1s.NT + g1s.NS -1.) / (g1s.NT * g1s.N1T) -
	   16.* (2.* g1s.NS + 2.* N +1.)/(g1s.NT * g1s.N1T) ) * (-0.5);
  PNSB = ( g1s.S1 * (536./9. + 8.* (2.* N + 1.) / (g1s.NS * g1s.N1S)) -
	   (16.* g1s.S1 + 52./3.- 8./(N * g1s.N1)) * g1s.S2 - 43./6. -
	   (151.* g1s.NFO + 263.* g1s.NT + 97.* g1s.NS + 3.* N + 9.) *
	   4./ (9.* g1s.NT * g1s.N1T) ) * (-0.5);
  PNSC = ( -160./9.* g1s.S1 + 32./3.* g1s.S2 + 4./3. +
	   16.*(11.*g1s.NS+5.*N-3.)/(9.* g1s.NS * g1s.N1S))*(-0.5);
}

void P1NSminus(dcomplex N, gamma1sums g1s, dcomplex &PNMA, dcomplex &PNSB, dcomplex &PNSC) {
  PNMA = ( 16.* g1s.S1 * (2.* N + 1.) / (g1s.NS * g1s.N1S) +
	   16.* (2.* g1s.S1 - 1./(N * g1s.N1)) * ( g1s.S2 - g1s.SSTR2M ) +
	   64.* g1s.SSCHLM + 24.* g1s.S2 - 3. - 8.* g1s.SSTR3M -
	   8.* (3.* g1s.NT + g1s.NS -1.) / (g1s.NT * g1s.N1T) +
	   16.* (2.* g1s.NS + 2.* N +1.) / (g1s.NT * g1s.N1T) ) * (-0.5); 
  PNSB = ( g1s.S1 * (536./9. + 8.* (2.* N + 1.) / (g1s.NS * g1s.N1S)) -
	   (16.* g1s.S1 + 52./3.- 8./(N * g1s.N1)) * g1s.S2 - 43./6. -
	   (151.* g1s.NFO + 263.* g1s.NT + 97.* g1s.NS + 3.* N + 9.) *
	   4./ (9.* g1s.NT * g1s.N1T) ) * (-0.5);
  PNSC = ( -160./9.* g1s.S1 + 32./3.* g1s.S2 + 4./3. +
	   16.*(11.*g1s.NS+5.*N-3.)/(9.* g1s.NS * g1s.N1S))*(-0.5);
}


//Singlet P1SG: PS (pure singlet) and QG...
void PS(dcomplex N, gamma1sums g1s, dcomplex &PPSA) {
  PPSA = (5.* g1s.NFI + 32.* g1s.NFO + 49.* g1s.NT+38.* g1s.NS + 28.* N + 8.) 
    / (g1s.NM * g1s.NT * g1s.N1T * g1s.N2S) * 2.;
}
   
void QG(dcomplex N, gamma1sums g1s, dcomplex &PQGA, dcomplex &PQGB) {
  PQGA = (-2.* g1s.S1 * g1s.S1 + 2.* g1s.S2 - 2.* g1s.SSTR2P) 
    * (g1s.NS + N + 2.) / (N * g1s.N1 * g1s.N2) 
    + (8.* g1s.S1 * (2.* N + 3.)) / (g1s.N1S * g1s.N2S)
    + 2.* (g1s.NN + 6.* g1s.NE + 15.* g1s.NSE + 25.* g1s.NSI + 36.* g1s.NFI
	   + 85.* g1s.NFO + 128.* g1s.NT + 104.* g1s.NS + 64.* N + 16.)
    / (g1s.NM * g1s.NT * g1s.N1T * g1s.N2T);
  PQGB = (2.* g1s.S1 * g1s.S1 - 2.* g1s.S2 + 5.) * (g1s.NS + N + 2.)
    / (N * g1s.N1 * g1s.N2) - 4.* g1s.S1 / g1s.NS
    + (11.* g1s.NFO + 26.* g1s.NT + 15.* g1s.NS + 8.* N + 4.)
    / (g1s.NT * g1s.N1T * g1s.N2) ;
}

void GQ(dcomplex N, gamma1sums g1s, dcomplex &PGQA, dcomplex &PGQB, dcomplex &PGQC) {
  PGQA = (- g1s.S1 * g1s.S1 + 5.* g1s.S1 - g1s.S2) * (g1s.NS + N + 2.) 
    / (g1s.NM * N * g1s.N1)  -  2.* g1s.S1 / g1s.N1S
    - (12.* g1s.NSI + 30.* g1s.NFI + 43.* g1s.NFO + 28.* g1s.NT - g1s.NS
       - 12.* N - 4.) / (2.* g1s.NM * g1s.NT * g1s.N1T) ;
  PGQB = (g1s.S1*g1s.S1 + g1s.S2 - g1s.SSTR2P) * (g1s.NS + N + 2.) / (g1s.NM * N * g1s.N1)
    - g1s.S1 * (17.* g1s.NFO + 41.* g1s.NS - 22.* N - 12.) 
    / (3.* g1s.NMS * g1s.NS * g1s.N1)
    + (109.* g1s.NN + 621.* g1s.NE + 1400.* g1s.NSE + 1678.* g1s.NSI
       + 695.* g1s.NFI - 1031.* g1s.NFO - 1304.* g1s.NT - 152.* g1s.NS
       + 432.* N + 144.) / (9.* g1s.NMS * g1s.NT * g1s.N1T * g1s.N2S);
  PGQC = (g1s.S1 - 8./3.) * (g1s.NS + N + 2.) / (g1s.NM * N * g1s.N1) + 1./ g1s.N1S;
  PGQC *= 4./3.;
}

void GG(dcomplex N, gamma1sums g1s, dcomplex &PGGA, dcomplex &PGGB, dcomplex &PGGC) {
  PGGA = - (2.* g1s.NFI + 5.* g1s.NFO + 8.* g1s.NT + 7.* g1s.NS- 2.* N - 2.)
    * 8.* g1s.S1 / (g1s.NMS * g1s.NS * g1s.N1S * g1s.N2S) -  67./9.* g1s.S1 + 8./3.
    - 4.* g1s.SSTR2P * (g1s.NS + N + 1.) / (g1s.NM * N * g1s.N1 * g1s.N2)
    + 2.* g1s.S1 * g1s.SSTR2P - 4.* g1s.SSCHLP + 0.5 * g1s.SSTR3P
    + (457.* g1s.NN + 2742.* g1s.NE + 6040.* g1s.NSE + 6098.* g1s.NSI
       + 1567.* g1s.NFI - 2344.* g1s.NFO - 1632.* g1s.NT + 560.* g1s.NS
       + 1488.* N + 576.) / (18.* g1s.NMS * g1s.NT * g1s.N1T * g1s.N2T);
  PGGB = (38.* g1s.NFO + 76.* g1s.NT + 94.* g1s.NS + 56.* N + 12.) *(-2.)
    / (9.* g1s.NM * g1s.NS * g1s.N1S * g1s.N2)  +  20./9.* g1s.S1  -  4./3.;
  PGGC= (2.*g1s.NSI + 4. * g1s.NFI + g1s.NFO - 10.* g1s.NT - 5. * g1s.NS - 4.* N
	 - 4.) * (-2.) / (g1s.NM * g1s.NT * g1s.N1T * g1s.N2) -1. ;
}



// ******** NON SINGLET ********
//
//Plus
dcomplex gamma1NSplus(dcomplex N, double nf, gamma1sums g1s){
  dcomplex PNPA, PNSB, PNSC;
  P1NSplus(N, g1s, PNPA, PNSB, PNSC);
  return CF *((CF-CA/2.)* PNPA + CA* PNSB + TF*(nf)* PNSC) /fourPi2;
}
dcomplex gamma1NSplus(dcomplex N, double nf){
  gamma1sums g1s;
  sums(N, g1s);
  return gamma1NSplus(N, nf, g1s);
}
//Minus=Valence
dcomplex gamma1NSminus(dcomplex N, double nf, gamma1sums g1s){
  dcomplex PNMA, PNSB, PNSC;
  P1NSminus(N, g1s, PNMA, PNSB, PNSC);
  return CF *((CF-CA/2.)* PNMA + CA* PNSB + TF*nf* PNSC) /fourPi2;
}
dcomplex gamma1NSminus(dcomplex N, double nf){
  gamma1sums g1s;
  sums(N, g1s);
  return gamma1NSminus(N, nf, g1s);
}


// ******** SINGLET *********
//
dcomplex gamma1SGqq(dcomplex N, double nf, gamma1sums g1s) {
  dcomplex PPSA;
  PS(N, g1s, PPSA);
  return gamma1NSplus(N, nf, g1s)  + TF*nf*CF*PPSA*4. /fourPi2;
}
dcomplex gamma1SGqg(dcomplex N, double nf, gamma1sums g1s) {
  dcomplex PQGA, PQGB;
  QG(N, g1s, PQGA, PQGB);
  return TF*nf * (CA * PQGA + CF * PQGB)*4. /fourPi2;
}
dcomplex gamma1SGgq(dcomplex N, double nf, gamma1sums g1s) {
  dcomplex PGQA, PGQB, PGQC;
  GQ(N, g1s, PGQA, PGQB, PGQC);
  return (CF*CF*PGQA + CF*CA*PGQB+TF*nf*CF*PGQC)*4. /fourPi2;
}
dcomplex gamma1SGgg(dcomplex N, double nf, gamma1sums g1s) {
  dcomplex PGGA, PGGB, PGGC;
  GG(N, g1s, PGGA, PGGB, PGGC);
  return (CA*CA*PGGA + TF*nf*(CA*PGGB+CF*PGGC))*4. /fourPi2;
}
//
dcomplex gamma1SGqq(dcomplex N, double nf) {
  gamma1sums g1s;
  sums(N, g1s);
  return gamma1SGqq(N, nf, g1s);
}
dcomplex gamma1SGqg(dcomplex N, double nf) {
  gamma1sums g1s;
  sums(N, g1s);
  return gamma1SGqg(N, nf, g1s);
}
dcomplex gamma1SGgq(dcomplex N, double nf) {
  gamma1sums g1s;
  sums(N, g1s);
  return gamma1SGgq(N, nf, g1s);
}
dcomplex gamma1SGgg(dcomplex N, double nf) {
  gamma1sums g1s;
  sums(N, g1s);
  return gamma1SGgg(N, nf, g1s);
}





// real versions
double gamma1SGqq(double N, double nf) {
  return real(gamma1SGqq(dcomplex(N,0.),nf));
}
double gamma1SGqg(double N, double nf) {
  return real(gamma1SGqg(dcomplex(N,0.),nf));
}
double gamma1SGgq(double N, double nf) {
  return real(gamma1SGgq(dcomplex(N,0.),nf));
}
double gamma1SGgg(double N, double nf) {
  return real(gamma1SGgg(dcomplex(N,0.),nf));
}
