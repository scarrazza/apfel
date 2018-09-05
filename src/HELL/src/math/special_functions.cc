#include <iostream>
#include <cmath>
#include "HELL/include/math/special_functions.hh"


using namespace std;

namespace HELLx {

//const double PISQ   =  9.86960440108935861883; 
const double PISQ6  =  1.64493406684822643647;

double ArcCsch(double a) { if(a<1e-4) return -log(a/2)+a*a/4; return log(1/a+sqrt(1+1/a/a)); }
double ArcCoth(double x) { return (log(1+1/x)-log(1-1/x))/2; }

double Li2(double x){
 double x_0 = -0.30;
 double x_1 = 0.25;
 double x_2 = 0.51;
 if (x == 1.) return PISQ6;
 if (x <= x_0){ 
   double temp = log(fabs(1.0-x));
   return -Li2(-x/(1.0-x)) - temp*temp/2 ; }
 else if (x < x_1){
   double z = - log(1.0-x);
   double temp = z*(1.0-z/4.0*(1.0-z/9.0*(1.0-z*z/100.0
                  *(1.0-5.0*z*z/294.0*(1.0-7.0*z*z/360.0
                  *(1.0-5.0*z*z/242.0*(1.0-7601.0*z*z/354900.0
                  *(1.0-91.0*z*z/4146.0*(1.0-3617.0*z*z/161840.0)
                   ))))))));
   return temp; }
   else if (x < x_2) return - Li2(-x) + Li2(x*x)/2.0 ;
   else { return PISQ6 - Li2(1.0-x) 
                  - log(fabs(x))*log(fabs(1.0-x)) ; }
}

double HPLmp(double x) {
  return Li2((1-x)/2) - Li2((1+x)/2) + (log(1-x)-log(1+x))*(log(1+x)+log(1-x)-log(4.))/2;
}


// the trilog, li3.
// good to ~ 10^(-14)
double Li3(double x){
 double x_0 = -1.0;
 double x_1 = -0.85;
 double x_2 = 0.25;
 double x_3 = 0.63;
 double x_4 =  1.0;
 if (x == 1.) return ZETA3;
 if (x == -1.) return - 0.75 * ZETA3;
 if (x <= x_0){ 
   double lnx = log(-x);
   return Li3(1.0/x) - PISQ6*lnx - lnx*lnx*lnx/6.0; }
 else if (x < x_1){
   return Li3(x*x)/4.0 - Li3(-x); }
   else if (x < x_2){
     double z = - log(1.0-x);
     double temp = z*(1.0-3.0*z/8.0*(1.0-17.0*z/81.0*(1.0-15*z/136.0
                    *(1.0-28.0*z/1875.0*(1.0+5.0*z/8.0*(1.0-304.0*z/7203.0
                    *(1.0+945.0*z/2432.0*(1.0-44.0*z/675.0*(1.0+7.0*z/24.0
                    *(1.0-26104.0*z/307461.0*(1.0+1925.0*z/8023.0
                    *(1.0-53598548.0*z/524808375.0
                    *(1.0+22232925.0*z/107197096.0
                     )))))))))))));
     return temp; }
     else if (x < x_3){
       return Li3(x*x)/4.0 - Li3(-x); }
       else if (x < x_4){
         double ln1x = log(1.0-x); 
         return -Li3(1.0-x) - Li3(-x/(1.0-x)) + ZETA3 + PISQ6*ln1x
	   - log(x)*ln1x*ln1x/2.0 + ln1x*ln1x*ln1x/6.0; }
       else { 
         double lnx = log(x);
         return Li3(1./x) + 2.0*PISQ6*lnx - lnx*lnx*lnx/6.0; }
}


// the tetralog, Li4.
// good to ~ 1.7 * 10^(-12) (worst areas:  x = 0.9501, also x = - 0.9701)
double Li4(double x){
 double x_0 = -1.0;
 double x_1 = -0.97;
 double x_2 = 0.25;
 double x_3 = 0.95;
 double x_4 =  1.0;
 if (x == -1) return -0.35 * PISQ6*PISQ6;
 if (x == 1) return 0.4 * PISQ6*PISQ6;
 if (x <= x_0){ 
   double lnx = log(-x);
   return - Li4(1./x) - 0.5 * PISQ6*lnx*lnx
          - 1./24. * lnx*lnx*lnx*lnx - 0.7 * PISQ6*PISQ6; }
 else if (x < x_1){
   return Li4(x*x)/8. - Li4(-x); }
   else if (x < x_2){
     double z = - log(1.-x);
     double temp = z*(1.-7.*z/16.*(1.-151.*z/567.*(1.0-411.*z/2416.
                    *(1.-24986.*z/256875.*(1.-805.*z/49972.
                    *(1.+583406.*z/1159683.*(1.-7455.*z/137272.
                    *(1.+659921.*z/2444175.*(1.-251559.*z/2639684.
                    *(1.+24259894.*z/136410197.*(1.-30625595.*z/218339046.
                    *(1.+2134239258113.*z/16698772722450.
                    *(1.-1640443805715.*z/8536957032452.
                     )))))))))))));
     return temp; }
     else if (x < x_3){
       return Li4(x*x)/8. - Li4(-x); }
       else if (x < x_4){
        double y = 1.-x; 
        double lny = log(y);
        return 0.4*PISQ6*PISQ6 + ZETA3 * log(1.-y) + 0.5*PISQ6 * y*y 
         + (-11./36.+0.5*PISQ6+1./6.*lny) * y*y*y
         + (11./24.*PISQ6+1./4.*lny-19./48.) * y*y*y*y
         + (7./24.*lny+5./12.*PISQ6-599./1440.) * y*y*y*y*y
         + (137./360.*PISQ6+5./16.*lny-79./192.) * y*y*y*y*y*y
	 + (7./20.*PISQ6-3343./8400.+29./90.*lny) * y*y*y*y*y*y*y
	  + (363./1120*PISQ6-21977./57600.+469./1440.*lny) * y*y*y*y*y*y*y*y;
                        }
       else { 
         double lnx = log(x);
         return - Li4(1/x) + PISQ6*lnx*lnx 
                - 1./24. * lnx*lnx*lnx*lnx + 0.8 * PISQ6*PISQ6; }
}


double dpsi(double Z, int M){
  int K1, K2;
  double SUB = 0. , SUBM;
  double ZZ = Z;
  if(true) { // if too close to the real axis...
  label1:
    SUBM = -1./ZZ;
    for(K1=1; K1<=M; K1++) {
      SUBM = - SUBM * K1 / ZZ;
    }
    if(ZZ<10.) { // ...use recurrence relation to push real(z) large enough
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
  double RZ = 1./ ZZ;
  double DZ = RZ * RZ;
  // SUB + asympt expansion (Abramowitz, Stengun, 6.4.11)
  return SUB + pow(-1.,M+1) * pow(RZ,M) * ( A1 + RZ * (A2 + RZ * (A3 + DZ * (A4 + DZ * (A5 + DZ * (A6 + A7 * DZ ))))) );
}

};

