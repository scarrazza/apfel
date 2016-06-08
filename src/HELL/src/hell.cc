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
    xx    = new double[Np1+Np2];
    xdPgg = new double[Np1+Np2];
    xdPqg = new double[Np1+Np2];
    for(int i=0; i<Np1+Np2; i++) {
      if(i<Np1) xx[i] = x_min * exp(i/(Np1-1.)*log(x_mid/x_min));
      else      xx[i] = x_mid + (i-Np1+1)*(x_max-x_mid)/(Np2-0.);
      infile >> xdPgg[i];
      if(isNLL) infile >> xdPqg[i];
    }
    infile.close();
  }

  double xTable::eval(double x) {
    if(x>1 || x<0) {
      cout << "\033[0;31m" << "Error: requesting resummed splitting function for unphysical value of x=" << x << " outside the physical range 0<x<=1" << "\033[0m" << endl;
      exit(45);
    }
    if(x>x_max) {
      cout << "\033[0;31m" << "Warning! Extrapolating out of interpolation range: x=" << x << " > x_max=" << x_max << "\033[0m" << endl;
      x = x_max;
    }
    if(x<x_min) {
      cout << "\033[0;31m" << "Warning! Extrapolating out of interpolation range: x=" << x << " < x_min=" << x_min << "\033[0m" << endl;
      x = x_min;
    }
    double ii;
    if(x<x_mid) ii = (Np1-1.)*log(x/x_min)/log(x_mid/x_min);
    else        ii = Np1-1.+Np2*(x-x_mid)/(x_max-x_mid);
    return ii;
  }
  void xTable::evalP(double x, double &dPgg, double &dPqg) {
    double ii = eval(x);
    int i = int(ii);
    if(i<0) {
      cout << "\033[0;31m" << "Error: this should never happen" << "\033[0m" << endl;
      abort();
    }
    //cout << setw(15) << x << setw(15) << xx[i] << setw(15) << xx[i+1] << setw(15) << xx[i] + (ii-i)*(xx[i+1]-xx[i]) << endl;
    dPgg   = xdPgg[i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdPgg[i+1]-xdPgg[i]));
    if(isNLL)
      dPqg = xdPqg[i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdPqg[i+1]-xdPqg[i]));
    return;
  }
  /*
  void xTable::evalC(double x, double &dC2g, double &dCLg) {
    double ii = eval(x);
    int i = int(ii);
    if(i<0) {
      cout << "\033[0;31m" << "Error: this should never happen" << "\033[0m" << endl;
      abort();
    }
    dC2g = xdC2g[i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdC2g[i+1]-xdC2g[i]));
    dCLg = xdCLg[i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdCLg[i+1]-xdCLg[i]));
    return;
  }
  */







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

  void HELLxnf::ReadTable(int k, map<int,xTable*> &T, string basename) {
    string sord = (_order == 1 ? "nlo" : "lo");
    ostringstream filename;
    itxT = T.find(k);
    if (itxT == T.end()) {
      filename << datapath_ << "/" << basename << "_" << sord << "_nf" << _nf << "_alphas" << sas(_alphas[k]) << ".table";
      T[k] = new xTable(filename.str(),_order);
      filename.str("");  filename.clear();
    }
    itxT = T.find(k+1);
    if (itxT == T.end()) {
      filename << datapath_ << "/" << basename << "_" << sord << "_nf" << _nf << "_alphas" << sas(_alphas[k+1]) << ".table";
      T[k+1] = new xTable(filename.str(),_order);
    }
  }


  sqmatrix<double> HELLxnf::DeltaP(double as, double x, Order matched_to_fixed_order) {
    double factor;
    int k = alphas_interpolation(as, _alphas, factor);
    ReadTable(k, xT, "xtable");
    double xdPgg[2], xdPqg[2];
    xT[k  ]->evalP(x,xdPgg[0],xdPqg[0]);
    xT[k+1]->evalP(x,xdPgg[1],xdPqg[1]);
    double dGgg = (xdPgg[0]+factor*(xdPgg[1]-xdPgg[0]))/x;
    double dGqg = (xdPqg[0]+factor*(xdPqg[1]-xdPqg[0]))/x;
    return sqmatrix<double>(dGgg, CF/CA*dGgg, dGqg, CF/CA*dGqg);
  }
  double HELLxnf::DeltaC(double as, double x, Order matched_to_fixed_order, string id) {
    if(_order == 0) return 0;
    double factor;
    int k = alphas_interpolation(as, _alphas, factor);
    ReadTable(k, xTC, "xtableC");
    double xdC2g[2], xdCLg[2];
    xTC[k  ]->evalP(x,xdC2g[0],xdCLg[0]);
    xTC[k+1]->evalP(x,xdC2g[1],xdCLg[1]);
    map<string,double> dCg;
    dCg["F2"] = (xdC2g[0]+factor*(xdC2g[1]-xdC2g[0]))/x;
    dCg["FL"] = (xdCLg[0]+factor*(xdCLg[1]-xdCLg[0]))/x;
    //return dCg[id];
    return dCg[id] / _nf;
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

  double HELLxnf::deltaMC2g (double as, double x, double m, Order matched_to_fixed_order) {
    //return DeltaC(as, x, matched_to_fixed_order, "F2m");
    return 0;
  }
  double HELLxnf::deltaMC2q (double as, double x, double m, Order matched_to_fixed_order) {
    return deltaMC2g(as, x, m, matched_to_fixed_order) * CF/CA;
  }

  double HELLxnf::deltaMCLg (double as, double x, double m, Order matched_to_fixed_order) {
    //return DeltaC(as, x, matched_to_fixed_order, "FLm");
    return 0;
  }
  double HELLxnf::deltaMCLq (double as, double x, double m, Order matched_to_fixed_order) {
    return deltaMCLg(as, x, m, matched_to_fixed_order) * CF/CA;
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

  // Delta P Matrix
  sqmatrix<double> HELLx::DeltaP(double as, double x,  Order matched_to_fixed_order) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->DeltaP(as, x, matched_to_fixed_order);
  }


  double HELLx::deltaC2g  (double as, double x, Order matched_to_fixed_order) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->deltaC2g(as, x, matched_to_fixed_order);
  }
  double HELLx::deltaC2q  (double as, double x, Order matched_to_fixed_order) {
    return deltaC2g(as, x, matched_to_fixed_order) * CF/CA;
  }

  double HELLx::deltaCLg  (double as, double x, Order matched_to_fixed_order) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->deltaCLg(as, x, matched_to_fixed_order);
  }
  double HELLx::deltaCLq  (double as, double x, Order matched_to_fixed_order) {
    return deltaCLg(as, x, matched_to_fixed_order) * CF/CA;
  }

  double HELLx::deltaMC2g (double as, double x, double m, Order matched_to_fixed_order) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->deltaMC2g(as, x, matched_to_fixed_order);
  }
  double HELLx::deltaMC2q (double as, double x, double m, Order matched_to_fixed_order) {
    return deltaMC2g(as, x, m, matched_to_fixed_order) * CF/CA;
  }

  double HELLx::deltaMCLg (double as, double x, double m, Order matched_to_fixed_order) {
    int nf = nf_of_as(as);
    return sxD[nf-nf_min]->deltaMCLg(as, x, matched_to_fixed_order);
  }
  double HELLx::deltaMCLq (double as, double x, double m, Order matched_to_fixed_order) {
    return deltaMCLg(as, x, m, matched_to_fixed_order) * CF/CA;
  }












};

