#ifndef __gammaNLO_hh__
#define __gammaNLO_hh__

#include "include/math/matrix.hh"


// LO
//------------------------------------------------------

dcomplex gamma0qq(dcomplex N);
dcomplex gamma0qg(dcomplex N, double nf);
dcomplex gamma0gq(dcomplex N);
dcomplex gamma0gg(dcomplex N, double nf);


// NLO
//------------------------------------------------------

struct gamma1sums {
  dcomplex S1, S2;
  dcomplex NS,NT,NFO,NFI,NSI,NSE,NE,NN;
  dcomplex N1,N2,NM,NMS,N1S,N1T,N2S,N2T;
  dcomplex N3,N4,N5,N6;
  dcomplex S11,S12,S13,S14,S15,S16;
  dcomplex SPMOM,SLC,SLV,SSCHLM,SSTR2M,SSTR3M,SSCHLP;
  dcomplex SSTR2P,SSTR3P;
};
extern void sums(dcomplex N, gamma1sums &g1s);


// NON SINGLET
//
//Plus
dcomplex gamma1NSplus(dcomplex N, double nf, gamma1sums g1s);
dcomplex gamma1NSplus(dcomplex N, double nf);
//Minus=Valence
dcomplex gamma1NSminus(dcomplex N, double nf, gamma1sums g1s);
dcomplex gamma1NSminus(dcomplex N, double nf);



//SINGLET
//
dcomplex gamma1SGqq(dcomplex N, double nf, gamma1sums g1s);
dcomplex gamma1SGqg(dcomplex N, double nf, gamma1sums g1s);
dcomplex gamma1SGgq(dcomplex N, double nf, gamma1sums g1s);
dcomplex gamma1SGgg(dcomplex N, double nf, gamma1sums g1s);
//
dcomplex gamma1SGqq(dcomplex N, double nf);
dcomplex gamma1SGqg(dcomplex N, double nf);
dcomplex gamma1SGgq(dcomplex N, double nf);
dcomplex gamma1SGgg(dcomplex N, double nf);


// real versions
double gamma1SGqq(double N, double nf);
double gamma1SGqg(double N, double nf);
double gamma1SGgq(double N, double nf);
double gamma1SGgg(double N, double nf);


#endif
