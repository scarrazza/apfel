/* -----------------------------------------

   _  _ ____ _    _    
   |__| |___ |    |    
   |  | |___ |___ |___ 
   
   HELL: High-Energy Large Logarithms

   Author: Marco Bonvini

   Computes Delta P_res = P_res - P_FixedOrder
   from a file containing parameters

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



namespace HELL {

  // **************************
  // ***   nf fixed part    ***
  // **************************

  HELLnf::HELLnf(int nf, LogOrder order, string prepath) {
    _order = order;
    string sord;
    if     (_order==0) {
      //sord = "/LO_nf";
      sord = "/NLO_nf";
      cout << "WARNING: LL resummation currently uses NLL resummed eigenvector (there are double counting problems with NLO or higher)" << endl;
    }
    else if(_order==1) sord = "/NLO_nf";
    else {
      cout << "ERROR: LogOrder can be only 0 or 1: NNLL or higher small-x resummation is NOT available." << endl;
      exit(0);
    }
    _nf = nf;
    //
    ostringstream os;
    os << nf;
    string snf = os.str();
    string filename = prepath+sord+string(snf)+".dat";
    ifstream infile(filename.c_str());
    if(!infile) {
      cout << "ERROR: No such file: " << filename << endl;
      cout << "       Maybe nf=" << nf << " is wrong?" << endl;
      exit(1);
    }
    if(_order==0)
      //infile >> Nas >> Np1 >> Np2;
      infile >> Nas >> Np1 >> Np2 >> Np3; /// !!!!! hack to be able to read NLL files as LL files!
    else
      infile >> Nas >> Np1 >> Np2 >> Np3;
    //
    alphas = new double[Nas];
    Npole1 = new double[Nas];
    Npole2 = new double[Nas];
    residue1 = new double[Nas];
    residue2 = new double[Nas];
    pDPp  = new double*[Nas];
    pDPqg = new double*[Nas];
    pDC2g = new double*[Nas];
    pdDC2g= new double*[Nas];
    //
    for(int i=0; i<Nas; i++) {
      pDPp  [i] = new double[Np1+Np2+1];
      pDPqg [i] = new double[Np3+2];
      pDC2g [i] = new double[Np3+2];
      pdDC2g[i] = new double[Np3+2];
      infile >> alphas[i] >> Npole1[i] >> Npole2[i] >> residue1[i] >> residue2[i];
      for(int j=0; j<=Np1+Np2; j++) infile >> pDPp[i][j];
      //if(_order>0) {
      if(_order>=0) { /// !!!!! hack to be able to read NLL files as LL files!
	for(int j=0; j<Np3+2; j++)    infile >> pDPqg[i][j];
	for(int j=0; j<Np3+2; j++)    infile >> pDC2g[i][j];
	for(int j=0; j<Np3+2; j++)    infile >> pdDC2g[i][j];
      }
    }
    infile.close();
  }

  HELLnf::~HELLnf() {
    delete alphas;
    delete Npole1;
    delete Npole2;
    delete residue1;
    delete residue2;
    delete[] pDPp;
    delete[] pDPqg;
    delete[] pDC2g;
    delete[] pdDC2g;
  }


  int HELLnf::GetOrder() {
    return _order;
  }
  int HELLnf::GetNf() {
    return _nf;
  }



  // factorial
  int fact[] = { 1, 1, 2, 6,
		 24,
		 120,
		 720,
		 5040,
		 40320,
		 362880 };


  // Delta P_+
  double HELLnf::xdeltaPplus_fit(double x, double *p) {
    double log1ox = -log(x);
    double res = 0.;
    for(int i=0; i<Np1; i++) {
      res += p[i] * pow(x, i+1);
    }
    for(int i=0; i<Np2; i++) {
      res += p[Np1+i] * pow(log1ox, i+1);
    }
    return res + p[Np1+Np2];
  }
  dcomplex HELLnf::deltaGammaplus_fit(dcomplex N, double *p) {
    dcomplex ooN = 1./N;
    dcomplex res = 0.;
    for(int i=0; i<Np1; i++) {
      res += p[i] / (N+1.+(double)i);
    }
    for(int i=0; i<Np2; i++) {
      res += p[Np1+i] * fact[i+1] * pow(ooN, i+2);
    }
    return res + p[Np1+Np2]*ooN;
  }
  // other functions
  double HELLnf::xdeltaFunc_fit_xspace(double x, double *p) {
    double log1ox = -log(x);
    double res = 0.;
    for(int i=0; i<Np3; i++) {
      res += p[i+2] * pow(log1ox, i);
    }
    return res + p[1]/pow(x,p[0]);;
  }
  dcomplex HELLnf::deltaFunc_fit_Nspace(dcomplex N, double *p) {
    dcomplex ooN = 1./N;
    dcomplex res = 0.;
    for(int i=0; i<Np3; i++) {
      res += p[i+2] * fact[i] * pow(ooN, i+1);
    }
    return res + p[1]/(N-p[0]);
  }

  int HELLnf::alphas_interpolation(double as, double &factor) {
    int k = 0;
    if(as < alphas[0] || as > alphas[Nas-1]) {
      cout << "ERROR: alpha_s=" << as << " out of interpolation range [" << alphas[0] << ", " << alphas[Nas-1] << "]" << endl;
      exit(22);
    }
    for(int i=1; i<Nas; i++) {
      if(as <= alphas[i]) break;
      k++;
    }
    //cout << alphas[k] << " < " << as << " < " << alphas[k+1] << endl;
    factor = (as-alphas[k]) / (alphas[k+1]-alphas[k]);
    return k;
  }
  int HELLnf::alphas_interpolation(double as, double &p1, double &p2, double &r1, double &r2, double &factor) {
    int k = alphas_interpolation(as, factor);
    p1 = Npole1[k] + factor * (Npole1[k+1] - Npole1[k]);
    p2 = Npole2[k] + factor * (Npole2[k+1] - Npole2[k]);
    r1 = residue1[k] + factor * (residue1[k+1] - residue1[k]);
    r2 = residue2[k] + factor * (residue2[k+1] - residue2[k]);
    return k;
  }



  // Delta P_+
  double HELLnf::xdeltaPplus(double as, double x) {
    double p1, p2, r1, r2, factor;
    int k = alphas_interpolation(as, p1, p2, r1, r2, factor);
    double dPf_k1 = xdeltaPplus_fit(x, pDPp[k]  );
    double dPf_k2 = xdeltaPplus_fit(x, pDPp[k+1]);
    return r1/pow(x, p1) + r2/pow(x, p2) + dPf_k1 + factor * (dPf_k2 - dPf_k1);
  }
  dcomplex HELLnf::deltaGammaplus(double as, dcomplex N, double *leadingPole) {
    N -= 1.;  // restores the common notation
    double p1, p2, r1, r2, factor;
    int k = alphas_interpolation(as, p1, p2, r1, r2, factor);
    dcomplex dgf_k1 = deltaGammaplus_fit(N, pDPp[k]  );
    dcomplex dgf_k2 = deltaGammaplus_fit(N, pDPp[k+1]);
    if(leadingPole != NULL) *leadingPole = p1+1;
    return r1/(N-p1) + r2/(N-p2) + dgf_k1 + factor * (dgf_k2 - dgf_k1);
  }
  double HELLnf::GammaplusNpole(double as) {
    double factor;
    int k = alphas_interpolation(as, factor);
    return 1 + Npole1[k] + factor * (Npole1[k+1] - Npole1[k]);
  }
  // Delta P_qg
  double HELLnf::xdeltaPqg(double as, double x) {
    if(_order==0) return 0;
    double factor;
    int k = alphas_interpolation(as, factor);
    double dPf_k1 = xdeltaFunc_fit_xspace(x, pDPqg[k]  );
    double dPf_k2 = xdeltaFunc_fit_xspace(x, pDPqg[k+1]);
    return dPf_k1 + factor * (dPf_k2 - dPf_k1);
  }
  dcomplex HELLnf::deltaGammaqg(double as, dcomplex N) {
    if(_order==0) return 0;
    N -= 1.;  // restores the common notation
    double factor;
    int k = alphas_interpolation(as, factor);
    dcomplex dgf_k1 = deltaFunc_fit_Nspace(N, pDPqg[k]  );
    dcomplex dgf_k2 = deltaFunc_fit_Nspace(N, pDPqg[k+1]);
    return dgf_k1 + factor * (dgf_k2 - dgf_k1);
  }
  // Delta C_2g
  double HELLnf::xdeltaC2g(double as, double x) {
    if(_order==0) return 0;
    double factor;
    int k = alphas_interpolation(as, factor);
    double dPf_k1 = xdeltaFunc_fit_xspace(x, pDC2g[k]  );
    double dPf_k2 = xdeltaFunc_fit_xspace(x, pDC2g[k+1]);
    return dPf_k1 + factor * (dPf_k2 - dPf_k1);
  }
  dcomplex HELLnf::deltaC2g_N(double as, dcomplex N) {
    if(_order==0) return 0;
    N -= 1.;  // restores the common notation
    double factor;
    int k = alphas_interpolation(as, factor);
    dcomplex dgf_k1 = deltaFunc_fit_Nspace(N, pDC2g[k]  );
    dcomplex dgf_k2 = deltaFunc_fit_Nspace(N, pDC2g[k+1]);
    return dgf_k1 + factor * (dgf_k2 - dgf_k1);
  }
  // Delta C_2g  ---  as-derivative
  double HELLnf::xdeltaC2g_deriv(double as, double x) {
    if(_order==0) return 0;
    double factor;
    int k = alphas_interpolation(as, factor);
    double dPf_k1 = xdeltaFunc_fit_xspace(x, pdDC2g[k]  );
    double dPf_k2 = xdeltaFunc_fit_xspace(x, pdDC2g[k+1]);
    return dPf_k1 + factor * (dPf_k2 - dPf_k1);
  }
  dcomplex HELLnf::deltaC2g_deriv_N(double as, dcomplex N) {
    if(_order==0) return 0;
    N -= 1.;  // restores the common notation
    double factor;
    int k = alphas_interpolation(as, factor);
    dcomplex dgf_k1 = deltaFunc_fit_Nspace(N, pdDC2g[k]  );
    dcomplex dgf_k2 = deltaFunc_fit_Nspace(N, pdDC2g[k+1]);
    return dgf_k1 + factor * (dgf_k2 - dgf_k1);
  }





  // Anomalous dimension Matrix
  sqmatrix<dcomplex> HELLnf::DeltaGammaLL(dcomplex N, double as, Order matched_to_fixed_order) {
    dcomplex dGgg = deltaGammaplus(as, N);
    if(matched_to_fixed_order < LO) dGgg += as*CA/M_PI/N;
    // Both of the following options are equally acceptable
    return sqmatrix<dcomplex>(dGgg, CF/CA*dGgg, 0., 0.);
    return sqmatrix<dcomplex>(dGgg, 0., 0., 0.);
  }
  sqmatrix<dcomplex> HELLnf::DeltaGammaNLL(dcomplex N, double as, Order matched_to_fixed_order) {
    dcomplex dGp  = deltaGammaplus(as, N);
    dcomplex dGqg = deltaGammaqg(as, N);
    dcomplex d0 = 0;
    if(matched_to_fixed_order < NLO) {
      dGp  += as*as/M_PI/M_PI/4. * _nf*(26.*CF-23.*CA)/9./(N-1.);
      dGqg += as*as/M_PI/M_PI/4. * _nf*20.*CA/9./(N-1.);
    }
    if(matched_to_fixed_order < LO)  {
      d0 = as/M_PI * _nf/3. * CF/CA;
      dGp  += as/M_PI * (CA/(N-1.) - (11.*CA+2.*_nf*(1.-2.*CF/CA))/12.) ;
      dGqg += as/M_PI * _nf/3.;
    }
    if(matched_to_fixed_order > NLO) {
      dGp  -= as*as*as/M_PI/M_PI/M_PI/8. * (-12.38818182*CA*CA*CA - 3.066013837*CA*CA*_nf + 6.132027674*CA*CF*_nf)/(N-1.)/(N-1.);
      dGqg -= as*as*as/M_PI/M_PI/M_PI/8. * 112.*CA*CA*_nf/27./(N-1.)/(N-1.);
    }
    dcomplex dGgg = dGp - CF/CA*dGqg;
    return sqmatrix<dcomplex>(dGgg, CF/CA*(dGgg+d0), dGqg, CF/CA*dGqg-d0);
  }
  sqmatrix<dcomplex> HELLnf::DeltaGamma(double as, dcomplex N, Order matched_to_fixed_order) {
    if(_order==0)
      return DeltaGammaLL(N, as, matched_to_fixed_order);
    else if(_order==1)
      return DeltaGammaNLL(N, as, matched_to_fixed_order);
    return sqmatrix<dcomplex>(0,0,0,0);
  }




  // Splitting Function Matrix
  sqmatrix<double> HELLnf::DeltaPLL(double x, double as, Order matched_to_fixed_order) {
    double dGgg = xdeltaPplus(as, x)/x;
    if(matched_to_fixed_order < LO) dGgg += as*CA/M_PI;
    // Both of the following options are equally acceptable
    return sqmatrix<double>(dGgg, CF/CA*dGgg, 0., 0.);
    return sqmatrix<double>(dGgg, 0., 0., 0.);
  }
  sqmatrix<double> HELLnf::DeltaPNLL(double x, double as, Order matched_to_fixed_order) {
    double dGp  = xdeltaPplus(as, x)/x;
    double dGqg = xdeltaPqg(as, x)/x;
    if(matched_to_fixed_order < NLO) {
      dGp  += as*as/M_PI/M_PI/4. * _nf*(26.*CF-23.*CA)/9./x;
      dGqg += as*as/M_PI/M_PI/4. * _nf*20.*CA/9./x;
    }
    if(matched_to_fixed_order < LO) {
      cout << "Error: cannot use pure NLL in x-space: delta(1-x) terms are missing" << endl;
      abort();
    }
    if(matched_to_fixed_order > NLO) {
      dGp  -= as*as*as/M_PI/M_PI/M_PI/8. * (-12.38818182*CA*CA*CA - 3.066013837*CA*CA*_nf + 6.132027674*CA*CF*_nf) * (-log(x)/x);
      dGqg -= as*as*as/M_PI/M_PI/M_PI/8. * 112.*CA*CA*_nf/27. * (-log(x)/x);
    }
    double dGgg = dGp - CF/CA*dGqg;
    return sqmatrix<double>(dGgg, CF/CA*dGgg, dGqg, CF/CA*dGqg);
  }
  sqmatrix<double> HELLnf::DeltaP(double as, double x, Order matched_to_fixed_order) {
    if(_order==0)
      return DeltaPLL(x, as, matched_to_fixed_order);
    else if(_order==1) {
      return DeltaPNLL(x, as, matched_to_fixed_order);
    }
    return sqmatrix<double>(0,0,0,0);
  }















  // **************************
  // ***  variable nf part  ***
  // **************************

  HELL::HELL(LogOrder order, string prepath) {
    _order = order;
    for(int inf=0; inf<values_of_nf; inf++)
      sxD[inf] = new HELLnf(nf_min+inf, order, prepath);
    //as_thr_init = false;
    init_mass_thresholds();
    init_as_thresholds();
  }
  HELL::~HELL() {
    for(int inf=0; inf<values_of_nf; inf++)
      delete sxD[inf];
  }

  //bool HELL::check_as_thr_init() {
  //  return as_thr_init;
  //}

  int HELL::GetOrder() {
    return _order;
  }


  // thresholds
  const double def_mass_c=1.21, def_mass_b=4.75, def_mass_t=173.; // default threshold masses
  const double mass_Z = 91.18, as_Z = 0.118;
  //
  double b0(int nf) {
    return (11.*CA-2.*nf)/12./M_PI;
  }
  double b1(int nf) {
    return (153.-19.*nf)/24./M_PI/M_PI/b0(nf);
  }
  double as_1loop(double mu, double alpha0, double m0, int nf) {
    return (alpha0/(1.+b0(nf)*alpha0*2.*log(mu/m0)));
  }
  double as_2loop(double mu, double alpha0, double m0, int nf) {
    return as_1loop(mu,alpha0,m0,nf)*(1-b1(nf)*as_1loop(mu,alpha0,m0,nf)*log(1.+b0(nf)*alpha0*2.*log(mu/m0)));
  }
  double HELL::as_LO(double mu) {
    double res=0;
    if(mu>mass_t) res=as_1loop(mu,as_LO(mass_t),mass_t,6);
    else if(mu<mass_b && mu>=mass_c) res=as_1loop(mu,as_LO(mass_b),mass_b,4);
    else if(mu<mass_c) res=as_1loop(mu,as_LO(mass_c),mass_c,3);
    else res=as_1loop(mu,as_Z,mass_Z,5);
    return res;
  }
  double HELL::as_NLO(double mu) {
    double res=0;
    if(mu>mass_t) res=as_2loop(mu,as_NLO(mass_t),mass_t,6);
    else if(mu<mass_b && mu>=mass_c) res=as_2loop(mu,as_NLO(mass_b),mass_b,4);
    else if(mu<mass_c) res=as_2loop(mu,as_NLO(mass_c),mass_c,3);
    else res=as_2loop(mu,as_Z,mass_Z,5);
    return res;
  }
  double HELL::as_thr(double mu) {
    if(_order==0) return as_LO(mu);
    else if(_order==1) return as_NLO(mu);
    return 0;
  }
  void HELL::init_mass_thresholds() {
    mass_c = def_mass_c;
    mass_b = def_mass_b;
    mass_t = def_mass_t;
  }
  void HELL::init_mass_thresholds(double mc, double mb, double mt) {
    mass_c = mc;
    mass_b = mb;
    mass_t = mt;
  }
  void HELL::init_as_thresholds() {
    as_c = as_thr(mass_c);
    as_b = as_thr(mass_b);
    as_t = as_thr(mass_t);
    //as_thr_init = true;
  }
  void HELL::init_as_thresholds(double as_of_mu(double)) {
    as_c = as_of_mu(mass_c);
    as_b = as_of_mu(mass_b);
    as_t = as_of_mu(mass_t);
    //as_thr_init = true;
  }
  void HELL::init_as_thresholds(double asc, double asb, double ast) {
    as_c = asc;
    as_b = asb;
    as_t = ast;
    //as_thr_init = true;
  }
  int HELL::nf_of_as(double as) {
    int res = 0;
    if(as<=as_t) res=6;
    else if(as>as_t && as<=as_b) res=5;
    else if(as>as_b && as<=as_c) res=4;
    else res=3;
    return res;
  }

  // Delta P_+
  double HELL::xdeltaPplus(double as, double x) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->xdeltaPplus(as,x);
  }
  dcomplex HELL::deltaGammaplus(double as, dcomplex N) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->deltaGammaplus(as, N);
  }
  // Delta P_qg
  double HELL::xdeltaPqg(double as, double x) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->xdeltaPqg(as, x);
  }
  dcomplex HELL::deltaGammaqg(double as, dcomplex N) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->deltaGammaqg(as, N);
  }
  // Delta C_2g
  double HELL::xdeltaC2g(double as, double x) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->xdeltaC2g(as, x);
  }
  dcomplex HELL::deltaC2g_N(double as, dcomplex N) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->deltaC2g_N(as, N);
  }
  // Delta C_2g  ---  as-derivative
  double HELL::xdeltaC2g_deriv(double as, double x) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->xdeltaC2g_deriv(as, x);
  }
  dcomplex HELL::deltaC2g_deriv_N(double as, dcomplex N) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->deltaC2g_deriv_N(as, N);
  }
  // Delta Gamma Matrix
  sqmatrix<dcomplex> HELL::DeltaGamma(double as, dcomplex N, Order matched_to_fixed_order) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->DeltaGamma(as, N, matched_to_fixed_order);
  }
  // Delta P Matrix
  sqmatrix<double> HELL::DeltaP(double as, double x,  Order matched_to_fixed_order) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->DeltaP(as, x, matched_to_fixed_order);
  }













};

