/* -----------------------------------------

   _  _ ____ _    _    
   |__| |___ |    |    
   |  | |___ |___ |___ x
   
   HELLx: High-Energy Large Logarithms - fast x-space version

   Author: Marco Bonvini

   Computes Delta P_res = P_res - P_FixedOrder
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


  xTable::xTable(string filename, bool nll) {
    isNLL = nll;
    ifstream infile(filename.c_str());
    if(!infile.good()) {
      cout << "Error reading table" << endl;
      abort();
    }
    infile >> Np1 >> Np2 >> x_min >> x_mid >> x_max;
    xx      = new double[Np1+Np2];
    xdPplus = new double[Np1+Np2];
    xdPqg   = new double[Np1+Np2];
    for(int i=0; i<Np1+Np2; i++) {
      if(i<Np1) xx[i] = x_min * exp(i/(Np1-1.)*log(x_mid/x_min));
      else      xx[i] = x_mid + (i-Np1+1)*(x_max-x_mid)/(Np2-0.);
      infile >> xdPplus[i];
      if(isNLL) infile >> xdPqg[i];
    }
    infile.close();
  }

  void xTable::eval(double x, double &dPplus, double &dPqg) {
    if(x>x_max) {
      x = x_max;
      cout << "warning: extrapolating out of interpolation range" << endl;
    }
    if(x<x_min) {
      x = x_min;
      cout << "warning: extrapolating out of interpolation range" << endl;
    }
    double ii;
    if(x<x_mid) ii = (Np1-1.)*log(x/x_min)/log(x_mid/x_min);
    else        ii = Np1-1.+Np2*(x-x_mid)/(x_max-x_mid);
    int i = int(ii);
    if(i<0) {
      cout << "error: this should never happen" << endl;
      abort();
    }
    //cout << setw(15) << x << setw(15) << xx[i] << setw(15) << xx[i+1] << setw(15) << xx[i] + (ii-i)*(xx[i+1]-xx[i]) << endl;
    dPplus = xdPplus[i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdPplus[i+1]-xdPplus[i]));
    if(isNLL)
      dPqg = xdPqg[i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdPqg[i+1]-xdPqg[i]));
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
    double astmp;
    while(info.good()) {
      info >> astmp;
      _alphas.push_back(astmp);
    }
    _alphas.pop_back();
  }

  HELLxnf::~HELLxnf() {
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
      cout << "ERROR: alpha_s=" << as << " out of interpolation range [" << vas[0] << ", " << vas[vas.size()-1] << "]" << endl;
      exit(22);
    }
    for(unsigned i=1; i<vas.size(); i++) {
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

  void HELLxnf::ReadTable(int k) {
    string sord = (_order == 1 ? "nlo" : "lo");
    ostringstream filename;
    itxT = xT.find(k);
    if (itxT == xT.end()) {
      filename << datapath_ << "/xtable_" << sord << "_nf" << _nf << "_alphas" << sas(_alphas[k]) << ".table";
      xT[k] = new xTable(filename.str(),_order);
      filename.str("");  filename.clear();
    }
    itxT = xT.find(k+1);
    if (itxT == xT.end()) {
      filename << datapath_ << "/xtable_" << sord << "_nf" << _nf << "_alphas" << sas(_alphas[k+1]) << ".table";
      xT[k+1] = new xTable(filename.str(),_order);
    }
  }




  // Delta P_+
  double HELLxnf::xdeltaPplus(double as, double x) {
    double factor;
    int k = alphas_interpolation(as, _alphas, factor);
    ReadTable(k);
    double xdPpl[2], xdPqg[2];
    xT[k  ]->eval(x,xdPpl[0],xdPqg[0]);
    xT[k+1]->eval(x,xdPpl[1],xdPqg[1]);
    double xdGp  = xdPpl[0]+factor*(xdPpl[1]-xdPpl[0]);
    return xdGp;
  }
  /*
  dcomplex HELLxnf::deltaGammaplus(double as, dcomplex N, double *leadingPole) {
    N -= 1.;  // restores the common notation
    double p1, p2, r1, r2, factor;
    int k = alphas_interpolation(as, p1, p2, r1, r2, factor);
    dcomplex dgf_k1 = deltaGammaplus_fit(N, pDPp[k]  );
    dcomplex dgf_k2 = deltaGammaplus_fit(N, pDPp[k+1]);
    if(leadingPole != NULL) *leadingPole = p1+1;
    return r1/(N-p1) + r2/(N-p2) + dgf_k1 + factor * (dgf_k2 - dgf_k1);
  }
  double HELLxnf::GammaplusNpole(double as) {
    double factor;
    int k = alphas_interpolation(as, factor);
    return 1 + Npole1[k] + factor * (Npole1[k+1] - Npole1[k]);
  }
  */
  // Delta P_qg
  double HELLxnf::xdeltaPqg(double as, double x) {
    if(_order==0) return 0;
    double factor;
    int k = alphas_interpolation(as, _alphas, factor);
    ReadTable(k);
    double xdPpl[2], xdPqg[2];
    xT[k  ]->eval(x,xdPpl[0],xdPqg[0]);
    xT[k+1]->eval(x,xdPpl[1],xdPqg[1]);
    double xdGqg = xdPqg[0]+factor*(xdPqg[1]-xdPqg[0]);
    return xdGqg;
  }
  /*
  dcomplex HELLxnf::deltaGammaqg(double as, dcomplex N) {
    if(_order==0) return 0;
    N -= 1.;  // restores the common notation
    double factor;
    int k = alphas_interpolation(as, factor);
    dcomplex dgf_k1 = deltaFunc_fit_Nspace(N, pDPqg[k]  );
    dcomplex dgf_k2 = deltaFunc_fit_Nspace(N, pDPqg[k+1]);
    return dgf_k1 + factor * (dgf_k2 - dgf_k1);
  }
  */




  /*
  // Anomalous dimension Matrix
  sqmatrix<dcomplex> HELLxnf::DeltaGammaLL(dcomplex N, double as, Order matched_to_fixed_order) {
    dcomplex dGgg = deltaGammaplus(as, N);
    if(matched_to_fixed_order < LO) dGgg += as*CA/M_PI/N;
    // Both of the following options are equally acceptable
    return sqmatrix<dcomplex>(dGgg, CF/CA*dGgg, 0., 0.);
    return sqmatrix<dcomplex>(dGgg, 0., 0., 0.);
  }
  sqmatrix<dcomplex> HELLxnf::DeltaGammaNLL(dcomplex N, double as, Order matched_to_fixed_order) {
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
  sqmatrix<dcomplex> HELLxnf::DeltaGamma(double as, dcomplex N, Order matched_to_fixed_order) {
    if(_order==0)
      return DeltaGammaLL(N, as, matched_to_fixed_order);
    else if(_order==1)
      return DeltaGammaNLL(N, as, matched_to_fixed_order);
    return sqmatrix<dcomplex>(0,0,0,0);
  }
  */



  // Splitting Function Matrix
  sqmatrix<double> HELLxnf::DeltaPLL(double x, double as, Order matched_to_fixed_order) {
    double dGgg = xdeltaPplus(as, x)/x;
    if(matched_to_fixed_order < LO) dGgg += as*CA/M_PI;
    // Both of the following options are equally acceptable
    return sqmatrix<double>(dGgg, CF/CA*dGgg, 0., 0.);
    return sqmatrix<double>(dGgg, 0., 0., 0.);
  }
  sqmatrix<double> HELLxnf::DeltaPNLL(double x, double as, Order matched_to_fixed_order) {
    double factor;
    int k = alphas_interpolation(as, _alphas, factor);
    ReadTable(k);
    double xdPpl[2], xdPqg[2];
    xT[k  ]->eval(x,xdPpl[0],xdPqg[0]);
    xT[k+1]->eval(x,xdPpl[1],xdPqg[1]);
    double dGp  = (xdPpl[0]+factor*(xdPpl[1]-xdPpl[0]))/x;
    double dGqg = (xdPqg[0]+factor*(xdPqg[1]-xdPqg[0]))/x;
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
  sqmatrix<double> HELLxnf::DeltaP(double as, double x, Order matched_to_fixed_order) {
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

  HELLx::HELLx(LogOrder order, string prepath) {
    _order = order;
    for(int inf=0; inf<values_of_nf; inf++)
      sxD[inf] = new HELLxnf(nf_min+inf, order, prepath);
    //as_thr_init = false;
    init_mass_thresholds();
    init_as_thresholds();
  }
  HELLx::~HELLx() {
    for(int inf=0; inf<values_of_nf; inf++)
      delete sxD[inf];
  }

  //bool HELLx::check_as_thr_init() {
  //  return as_thr_init;
  //}

  int HELLx::GetOrder() {
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
  double HELLx::as_LO(double mu) {
    double res=0;
    if(mu>mass_t) res=as_1loop(mu,as_LO(mass_t),mass_t,6);
    else if(mu<mass_b && mu>=mass_c) res=as_1loop(mu,as_LO(mass_b),mass_b,4);
    else if(mu<mass_c) res=as_1loop(mu,as_LO(mass_c),mass_c,3);
    else res=as_1loop(mu,as_Z,mass_Z,5);
    return res;
  }
  double HELLx::as_NLO(double mu) {
    double res=0;
    if(mu>mass_t) res=as_2loop(mu,as_NLO(mass_t),mass_t,6);
    else if(mu<mass_b && mu>=mass_c) res=as_2loop(mu,as_NLO(mass_b),mass_b,4);
    else if(mu<mass_c) res=as_2loop(mu,as_NLO(mass_c),mass_c,3);
    else res=as_2loop(mu,as_Z,mass_Z,5);
    return res;
  }
  double HELLx::as_thr(double mu) {
    if(_order==0) return as_LO(mu);
    else if(_order==1) return as_NLO(mu);
    return 0;
  }
  void HELLx::init_mass_thresholds() {
    mass_c = def_mass_c;
    mass_b = def_mass_b;
    mass_t = def_mass_t;
  }
  void HELLx::init_mass_thresholds(double mc, double mb, double mt) {
    mass_c = mc;
    mass_b = mb;
    mass_t = mt;
  }
  void HELLx::init_as_thresholds() {
    as_c = as_thr(mass_c);
    as_b = as_thr(mass_b);
    as_t = as_thr(mass_t);
    //as_thr_init = true;
  }
  void HELLx::init_as_thresholds(double as_of_mu(double)) {
    as_c = as_of_mu(mass_c);
    as_b = as_of_mu(mass_b);
    as_t = as_of_mu(mass_t);
    //as_thr_init = true;
  }
  void HELLx::init_as_thresholds(double asc, double asb, double ast) {
    as_c = asc;
    as_b = asb;
    as_t = ast;
    //as_thr_init = true;
  }
  int HELLx::nf_of_as(double as) {
    int res = 0;
    if(as<=as_t) res=6;
    else if(as>as_t && as<=as_b) res=5;
    else if(as>as_b && as<=as_c) res=4;
    else res=3;
    return res;
  }

  // Delta P_+
  double HELLx::xdeltaPplus(double as, double x) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->xdeltaPplus(as,x);
  }
  //dcomplex HELLx::deltaGammaplus(double as, dcomplex N) {
  //  int nf = nf_of_as(as);
  //  return sxD[nf-nf_min]->deltaGammaplus(as, N);
  //}
  // Delta P_qg
  double HELLx::xdeltaPqg(double as, double x) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->xdeltaPqg(as, x);
  }
  //dcomplex HELLx::deltaGammaqg(double as, dcomplex N) {
  //  int nf = nf_of_as(as);
  //  return sxD[nf-nf_min]->deltaGammaqg(as, N);
  //}
  /*
  // Delta C_2g
  double HELLx::xdeltaC2g(double as, double x) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->xdeltaC2g(as, x);
  }
  dcomplex HELLx::deltaC2g_N(double as, dcomplex N) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->deltaC2g_N(as, N);
  }
  // Delta C_2g  ---  as-derivative
  double HELLx::xdeltaC2g_deriv(double as, double x) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->xdeltaC2g_deriv(as, x);
  }
  dcomplex HELLx::deltaC2g_deriv_N(double as, dcomplex N) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->deltaC2g_deriv_N(as, N);
  }
  */
  // Delta P Matrix
  sqmatrix<double> HELLx::DeltaP(double as, double x,  Order matched_to_fixed_order) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->DeltaP(as, x, matched_to_fixed_order);
  }
  // Delta Gamma Matrix
  //sqmatrix<dcomplex> HELLx::DeltaGamma(double as, dcomplex N, Order matched_to_fixed_order) {
  //  int nf = nf_of_as(as);
  //  return sxD[nf-nf_min]->DeltaGamma(as, N, matched_to_fixed_order);
  //}













};

