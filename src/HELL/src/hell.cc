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

#include "HELL/include/hell.hh"


using namespace std;



namespace HELLx {

  const double ZETA2 = 1.6449340668482264;
  const double ZETA3 = 1.2020569031595942855;
  const double ZETA4 = 1.082323233711138191516;
  double Li2(double x){
    double x_0 = -0.30;
    double x_1 = 0.25;
    double x_2 = 0.51;
    if (x == 1.) return ZETA2;
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
    else { return ZETA2 - Li2(1.0-x) 
	- log(fabs(x))*log(fabs(1.0-x)) ; }
  }
  double HPLmp(double x) {
    return Li2((1-x)/2) - Li2((1+x)/2) + (log(1-x)-log(1+x))*(log(1+x)+log(1-x)-log(4.))/2;
  }
  inline double ArcCsch(double a) { if(a<1e-4) return -log(a/2)+a*a/4; return log(1/a+sqrt(1+1/a/a)); }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int _damping = 2, _dampingsqrt = 4;
  const double c0 = CA/M_PI*4*log(2.);
  const double k0 = CA/M_PI*28.*ZETA3;
  double a11()       { return CA/M_PI; }
  double a10(int nf) { return -(11.*CA + 2.*nf*(1.-2.*CF/CA))/12./M_PI; }
  double a21(int nf) { return nf*(26*CF-23*CA)/36/M_PI/M_PI; }
  double beta0(int nf) { return (11*CA-2*nf)/12./M_PI; }
  double fmom(double N) { return 4*N/(N+1)/(N+1); }
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
  double d1(int k, int j) {
    double res = 0;
    for(int i=0; i<=k; i++) {
      res += binomial(k,i) * pow(-1.,i) /(2.+i) /binomial(4+2*i+j,j);
    }
    return res;
  }
  double gammagg2(double N, int nf) {
    return beta0(nf)* (3*k0/32.-c0) * (1./N - fmom(N));
  }
  double gammagg3(double N, int nf) {
    double res = 0;
    double a1 = a11(), a0=a10(nf);
    double A = ( CA*CA*(-74./27. +11./12.*ZETA2 +5./2.*ZETA3) +nf*((4*CA+7*CF)/27. +(CA-2*CF)*ZETA2/6.) )/M_PI/M_PI;
    double chi11 = a10(nf)*a11();
    double chi02 = a21(nf);
    double gamma0p = -a1/N/N - 2*(a1+a0)/(N+1)/(N+1);
    res  = -gamma0p * ( A + chi11/(N+1.) + chi02/pow(N+1.,2) - (A+chi11/2.+chi02/4.)*fmom(N) );
    res += CA/M_PI* ( a1*a1*  ( 3*ZETA4/N -2*ZETA3/N/N -3/(N+1) -1/(N+1)/(N+1) )
		      +a1*a0* (-2*ZETA3/N +2/(N+1) +1/(N+1)/(N+1) )
		      +4*pow(a0 + a1,2)*              ( 3*ZETA4/(N+1) -2*ZETA3/(N+1)/(N+1) +0*1/pow(N+1,5) )
		      -2*(a0*a0 + 4*a0*a1 + 3*a1*a1)* ( -2*ZETA3/(N+1) +1/pow(N+1,4) )
		      );
    double b0 = beta0(nf);
    res += b0*b0*k0/16.* (1/N - fmom(N));
    double gamma0LL = a1/N + a0/(N+1);
    double gamma1LL = b0*(3*k0/32.-c0) * (1/N - fmom(N));
    double h1 = 5./3., h2 = 14./9.;
    res += (1-CF/CA) * (h2*gamma0LL*(gamma0LL-b0) + h1*gamma1LL) * nf/3./M_PI;
    return res;
  }
  double Ppl3(double x, int nf, bool domomcons=true) {
    double res = 0;
    double b0 = beta0(nf);
    double A = ( CA*CA*(-74./27. +11./12.*ZETA2 +5./2.*ZETA3) +nf*((4*CA+7*CF)/27. +(CA-2*CF)*ZETA2/6.) )/M_PI/M_PI;
    double a1 = a11(), a0=a10(nf);
    double chi11 = a10(nf)*a11();
    double chi02 = a21(nf);
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
    res += b0*b0*k0/16* (1/x-4*(1+log(x)));
    // momentum conservation
    double momcons = 0;
    int k=_damping, j=_dampingsqrt;
    if(_damping==2 && _dampingsqrt==4) {
      double coeffs[9] = {1, -4, 4, 4, -10, 4, 4, -4, 1};
      for(int n=0; n<=8; n++) {
	  momcons += coeffs[n] * gammagg3(1.+n/2.,nf);
      }
    } else {
      for(int n=0; n<=k; n++) {
	for(int m=0; m<=j; m++) {
	  momcons += binomial(k,n) * binomial(j,m) * pow(-1.,n+m) * gammagg3(1.+n+m/2.,nf);
	}
      }
    }
    momcons /= d1(k,j);
    if(domomcons) res -= momcons;
    return res * pow(1-x,k)*pow(1-sqrt(x),j);
  }
  double Pqg3(double x, int nf) {
    double res = 0;
    double h1 = 5./3., h2 = 14./9.;
    double b0 = beta0(nf);
    double A = b0*(3*k0/32-c0);
    //double A = b0*(-c0/4);
    double a1 = a11(), a0=a10(nf);
    res = h2*a1*a1*Pole(x,2,0) + (h2*a0*a0+4*h1*A)*Pole(x,0,2) + 2*h2*a1*a0*Pole(x,1,1) + (h1*A-b0*a1*h2)*Pole(x,1,0) - (h1*4*A+b0*a0)*Pole(x,0,1);
    int k=_damping, j=_dampingsqrt;
    return res*nf/3./M_PI * pow(1-x,k)*pow(1-sqrt(x),j);
  }
  double Pgg2(double x, int nf, bool domomcons=true) {
    double res = beta0(nf)*(3*k0/32.-c0)* (1/x-4*(1+log(x)));
    // momentum conservation
    double momcons = 0;
    int k=_damping, j=_dampingsqrt;
    for(int n=0; n<=k; n++) {
      for(int m=0; m<=j; m++) {
	momcons += binomial(k,n) * binomial(j,m) * pow(-1.,n+m) * gammagg2(1.+n+m/2.,nf);
      }
    }
    momcons /= d1(k,j);
    if(domomcons) res -= momcons;
    return res * pow(1-x,k)*pow(1-sqrt(x),j);
  }
  double Pplus0(double x, int nf) {
    double a1 = a11(), a0=a10(nf);
    double res = a1/x + a0;
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k) * pow(1-sqrt(x),j);  // here j was 6, now it is 4
  }
  double Pplus0sq(double x, int nf) {
    double a1 = a11(), a0=a10(nf);
    double res = a1*a1*Pole(x,2,0) + a0*a0*Pole(x,0,2) + 2*a1*a0*Pole(x,1,1);
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k) * pow(1-sqrt(x),j);  // here j was 6, now it is 4
  }
  double Pplus1(double x, int nf) {
    double res = beta0(nf)*(3*k0/32.-c0) * (1/x-4*(1+log(x)));
    int k=_damping, j=_dampingsqrt;
    return res * pow(1-x,k) * pow(1-sqrt(x),j);
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
      //cout << "\033[0;31m" << "HELLx: Warning! Extrapolating out of interpolation range: m/Q=" << mQ << " > " << mQvals[Nmass-1] << " for as=" << as << ", nf=" << nf << "\033[0m" << endl;
    }
    if(mQ<mQvals[0]) {
      //cout << "\033[0;31m" << "HELLx: Warning! Extrapolating out of interpolation range: m/Q=" << mQ << " < " << mQvals[0]       << " for as=" << as << ", nf=" << nf << "\033[0m" << endl;
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
    mQvals = new double [Nmass];
    xdKhg  = new double*[Nmass];
    xdC2g  = new double*[Nmass];
    xdCLg  = new double*[Nmass];
    for(int m=0; m<Nmass; m++) {
      xdKhg[m] = new double[Np1+Np2];
      xdC2g[m] = new double[Np1+Np2];
      xdCLg[m] = new double[Np1+Np2];
      *infile >> mQvals[m];
      for(int i=0; i<Np1+Np2; i++) {
	*infile >> xdKhg[m][i] >> xdC2g[m][i] >> xdCLg[m][i];
      }
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
  void xTableCm::eval(double x, double mQ, double &dKhg, double &dC2g, double &dCLg, double as, int nf) {
    double ii = interpolate(x);
    int i = int(ii);
    if(i<0 || ii<i) {
      cout << "\033[0;31m" << "HELLx: Error! This should never happen" << "\033[0m" << endl;
      abort();
    }
    double *mdKhg, *mdC2g, *mdCLg;
    mdKhg = new double[Nmass];
    mdC2g = new double[Nmass];
    mdCLg = new double[Nmass];
    for(int m=0; m<Nmass; m++) {
      mdKhg[m] = xdKhg[m][i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdKhg[m][i+1]-xdKhg[m][i]));
      mdC2g[m] = xdC2g[m][i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdC2g[m][i+1]-xdC2g[m][i]));
      mdCLg[m] = xdCLg[m][i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdCLg[m][i+1]-xdCLg[m][i]));
    }
    dKhg = minterpolate(mQ, mQvals, mdKhg, Nmass, x, as, nf);
    dC2g = minterpolate(mQ, mQvals, mdC2g, Nmass, x, as, nf);
    dCLg = minterpolate(mQ, mQvals, mdCLg, Nmass, x, as, nf);
    delete[] mdKhg;
    delete[] mdC2g;
    delete[] mdCLg;
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
    filename << datapath << "/lo_nf" << _nf << ".info";
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
    deleteTable(xT);
    deleteTable(xTC);
    deleteTable(xTCm);
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
      if(as <= vas[i]) break;
      k++;
    }
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
    for(int i=k; i<k+2; i++) {
      itxT = T.find(i);
      if (itxT == T.end()) {
	filename.str("");  filename.clear();
	filename << datapath_ << "/" << basename << "_" << sord << "_nf" << _nf << "_alphas" << sas(_alphas[i]) << ".table";
	T[i] = new S(filename.str(), _order);
      }
    }
  }


  sqmatrix<double> HELLxnf::DeltaP(double as, double x, Order matched_to_fixed_order) {
    if(_order==1 && (matched_to_fixed_order < NLO || matched_to_fixed_order > NNLO)) {
      cout << "Error: trying to match NLL resummation in splitting functions to a fixed order different from NLO or NNLO" << endl;
      exit(45);
    }
    if(_order==0 && (matched_to_fixed_order < LO || matched_to_fixed_order > NLO)) {
      cout << "Error: trying to match LL resummation in splitting functions to a fixed order different from LO or NLO" << endl;
      exit(45);
    }
    double factor;
    int k = alphas_interpolation(as, _alphas, factor);
    ReadTable(k, xT, "xtable");
    double xdPgg[2], xdPqg[2];
    xT[k  ]->eval(x,xdPgg[0],xdPqg[0]);
    xT[k+1]->eval(x,xdPgg[1],xdPqg[1]);
    double dGgg = (xdPgg[0]+factor*(xdPgg[1]-xdPgg[0]))/x;
    double dGqg = (xdPqg[0]+factor*(xdPqg[1]-xdPqg[0]))/x;
    if(_order==0 && matched_to_fixed_order == LO) {
      dGgg += as*as * Pgg2(x,_nf);
    }
    if(_order==1 && matched_to_fixed_order == NLO) {
      dGgg += as*as*as * ( Ppl3(x,_nf) - CF/CA*Pqg3(x,_nf) );
      dGqg += as*as*as *   Pqg3(x,_nf);
    }
    return sqmatrix<double>(dGgg, CF/CA*dGgg, dGqg, CF/CA*dGqg);
  }
  double HELLxnf::DeltaC(double as, double x, Order matched_to_fixed_order, string id) {
    if(_order==0) return 0;
    if(matched_to_fixed_order < NLO || matched_to_fixed_order > N3LO) {
      cout << "Error: trying to match resummed massless coefficient functions to a fixed order different from NLO or NNLO" << endl;
      exit(46);
    }
    double factor;
    int k = alphas_interpolation(as, _alphas, factor);
    ReadTable(k, xTC, "xtableC");
    //ReadTable(k, xTCm, "xtableCm");
    double xdC2g[2], xdCLg[2];
    xTC[k  ]->eval(x,xdC2g[0],xdCLg[0]);
    xTC[k+1]->eval(x,xdC2g[1],xdCLg[1]);
    map<string,double> dCg;
    dCg["F2"] = (xdC2g[0]+factor*(xdC2g[1]-xdC2g[0]))/x;
    dCg["FL"] = (xdCLg[0]+factor*(xdCLg[1]-xdCLg[0]))/x;
    if(matched_to_fixed_order == NLO) {
      double P0 = Pplus0(x,_nf);
      double h12 = 43./18.-ZETA2;
      double h1L = -1./3;
      dCg["F2"] += as*as * P0 * h12*_nf/3./M_PI;
      dCg["FL"] += as*as * P0 * h1L*_nf/3./M_PI;
    }
    if(matched_to_fixed_order == N3LO) {
      double P0sq = Pplus0sq(x,_nf);
      double P0   = Pplus0  (x,_nf);
      double P1   = Pplus1  (x,_nf);
      double h12 = 43./18.-ZETA2;
      double h1L = -1./3;
      double h22 = 3.2982926655872 /2.;
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
    ReadTable(k, xTCm, "xtableCm");
    double xdKhg[2], xdC2g[2], xdCLg[2];
    // converts m/Q-interpolation into m-interpolation
    double mQ0 = mQ * Qofalphas(as,_nf)/Qofalphas(_alphas[k  ],_nf);
    double mQ1 = mQ * Qofalphas(as,_nf)/Qofalphas(_alphas[k+1],_nf);
    xTCm[k  ]->eval(x, mQ0, xdKhg[0],xdC2g[0],xdCLg[0],as,_nf);
    xTCm[k+1]->eval(x, mQ1, xdKhg[1],xdC2g[1],xdCLg[1],as,_nf);
    //cout << setw(14) << mQ
    //	 << setw(14) << xdKhg[0] << setw(14) << xdC2g[0] << setw(14) << xdCLg[0]
    //	 << setw(14) << xdKhg[1] << setw(14) << xdC2g[1] << setw(14) << xdCLg[1] << endl;
    map<string,double> dCg;
    dCg["Khg"] = (xdKhg[0]+factor*(xdKhg[1]-xdKhg[0]))/x;
    dCg["F2m"] = (xdC2g[0]+factor*(xdC2g[1]-xdC2g[0]))/x;
    dCg["FLm"] = (xdCLg[0]+factor*(xdCLg[1]-xdCLg[0]))/x;
    if(matched_to_fixed_order == NLO) {
      double P0 = Pplus0(x,_nf);
      double mQ2 = mQ*mQ;
      double lmQ2 = log(mQ2);
      double sq = sqrt(1+4*mQ2);
      double h1K = -(28.+30.*lmQ2+9.*lmQ2*lmQ2)/18.;
      double h12 = (5+3*lmQ2)/6 - (1-mQ2)*HPLmp(-1./sq)/sq + ArcCsch(2*sqrt(mQ2))*(13-10*mQ2+6*(1-mQ2)*lmQ2)/3/sq;
      double h1L = (-1+12*mQ2+3*(1+6*mQ2)*lmQ2)/3/(1+4*mQ2) + 4*mQ2*(1+3*mQ2)*HPLmp(-1./sq)/(1+4*mQ2)/sq + ArcCsch(2*sqrt(mQ2))*(6+8*mQ2*(1-6*mQ2)-24*mQ2*(1+3*mQ2)*lmQ2)/3/(1+4*mQ2)/sq;
      dCg["Khg"] += as*as * P0 * h1K/3./M_PI;
      dCg["F2m"] += as*as * P0 * h12/3./M_PI;
      dCg["FLm"] += as*as * P0 * h1L/3./M_PI;
    }
    return dCg[id];
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
    if(4*m_Q_ratio*x/(1-x)>1) return 0;
    if(m_Q_ratio < 0.002) return deltaC2g(as, x, matched_to_fixed_order)/_nf + deltaKhg(as, x, m_Q_ratio, matched_to_fixed_order);
    return DeltaCm(as, x, matched_to_fixed_order, "F2m", m_Q_ratio);
  }
  double HELLxnf::deltaMC2q (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaMC2g(as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }

  double HELLxnf::deltaMCLg (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    if(4*m_Q_ratio*x/(1-x)>1) return 0;
    if(m_Q_ratio < 0.002) return deltaCLg(as, x, matched_to_fixed_order)/_nf;
    return DeltaCm(as, x, matched_to_fixed_order, "FLm", m_Q_ratio);
  }
  double HELLxnf::deltaMCLq (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaMCLg(as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
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
      cout << "HELLx: You requested matching function nf=6 scheme. Isn't it too much?" << endl;
      return 0;//deltaC2g(nf, as, x, matched_to_fixed_order);
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












};

