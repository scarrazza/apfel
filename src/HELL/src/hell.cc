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

  xTable::xTable(string filename, bool nll) {
    isNLL = nll;
    ifstream infile(filename.c_str());
    if(!infile.good()) {
      cout << "\033[0;31m" << "HELLx: Error reading table" << "\033[0m" << endl;
      exit(0);
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
      cout << "\033[0;31m" << "HELLx: Error: requesting resummed splitting function for unphysical value of x=" << x << " outside the physical range 0<x<=1" << "\033[0m" << endl;
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
  void xTable::evalP(double x, double &dPgg, double &dPqg) {
    double ii = eval(x);
    int i = int(ii);
    if(i<0) {
      cout << "\033[0;31m" << "HELLx: Error: this should never happen" << "\033[0m" << endl;
      abort();
    }
    dPgg = xdPgg[i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdPgg[i+1]-xdPgg[i]));
    dPqg = 0;
    if(isNLL)
      dPqg = xdPqg[i] + (i==Np1+Np2-1 ? 0 : (ii-i)*(xdPqg[i+1]-xdPqg[i]));
    return;
  }
  /*
  void xTable::evalC(double x, double &dC2g, double &dCLg) {
    double ii = eval(x);
    int i = int(ii);
    if(i<0) {
      cout << "\033[0;31m" << "HELLx: Error: this should never happen" << "\033[0m" << endl;
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
      cout << "\033[0;31m" << "HELLx: ERROR: alpha_s=" << as << " out of interpolation range [" << vas[0] << ", " << vas[vas.size()-1] << "]" << "\033[0m" << endl;
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
    if(_order==1 && matched_to_fixed_order != NLO) {
      cout << "Error: trying to match NLL resummation in splitting functions to a fixed order different from NLO" << endl;
      exit(45);
    }
    if(_order==0 && matched_to_fixed_order != LO) {
      cout << "Error: trying to match LL resummation in splitting functions to a fixed order different from LO" << endl;
      exit(45);
    }
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
    if(_order==0) return 0;
    double factor;
    int k = alphas_interpolation(as, _alphas, factor);
    ReadTable(k, xTC, "xtableC");
    //ReadTable(k, xTCm, "xtableCm");
    double xdC2g[2], xdCLg[2];
    xTC[k  ]->evalP(x,xdC2g[0],xdCLg[0]);
    xTC[k+1]->evalP(x,xdC2g[1],xdCLg[1]);
    map<string,double> dCg;
    dCg["F2"] = (xdC2g[0]+factor*(xdC2g[1]-xdC2g[0]))/x;
    dCg["FL"] = (xdCLg[0]+factor*(xdCLg[1]-xdCLg[0]))/x;
    if(matched_to_fixed_order == NLO) {
      double fixedorderpole = as/M_PI * (CA/x - (11*CA+2*_nf*(1-2*CF/CA))/12.) *as/M_PI *_nf/3.  * pow(1-x,2.) * pow(1-sqrt(x),6.);
      //double fixedorderpole = as/M_PI *  CA/x                                  *as/M_PI *_nf/3.  * pow(1-x,2.) * pow(1-sqrt(x),6.);
      dCg["F2"] += fixedorderpole * (43./18.-ZETA2);
      dCg["FL"] += fixedorderpole * (-1./3.);
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

  double HELLxnf::deltaMC2g (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    //return DeltaC(as, x, matched_to_fixed_order, "F2m");
    return 0;
  }
  double HELLxnf::deltaMC2q (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    return deltaMC2g(as, x, m_Q_ratio, matched_to_fixed_order) * CF/CA;
  }

  double HELLxnf::deltaMCLg (double as, double x, double m_Q_ratio, Order matched_to_fixed_order) {
    //return DeltaC(as, x, matched_to_fixed_order, "FLm");
    return 0;
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
      cout << "\033[0;31m" << "HELLx: Error: nf out of range " << nf << "\033[0m" << endl;
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
    check_nf(nf);
    return sxD[nf-3]->deltaC2g(as, x, matched_to_fixed_order) * CF/CA;
  }
  //double HELLx::deltaC2q  (int nf, double as, double x, Order matched_to_fixed_order) {
  //  return deltaC2g(as, x, matched_to_fixed_order) * CF/CA;
  //}

  double HELLx::deltaCLg  (int nf, double as, double x, Order matched_to_fixed_order) {
    check_nf(nf);
    return sxD[nf-3]->deltaCLg(as, x, matched_to_fixed_order);
  }
  double HELLx::deltaCLq  (int nf, double as, double x, Order matched_to_fixed_order) {
    check_nf(nf);
    return sxD[nf-3]->deltaCLg(as, x, matched_to_fixed_order) * CF/CA;
  }
  //double HELLx::deltaCLq  (int nf, double as, double x, Order matched_to_fixed_order) {
  //  return deltaCLg(as, x, matched_to_fixed_order) * CF/CA;
  //}

  double HELLx::deltaMC2g (int nf, double as, double x, double m, Order matched_to_fixed_order) {
    check_nf(nf);
    return sxD[nf-3]->deltaMC2g(as, x, matched_to_fixed_order);
  }
  double HELLx::deltaMC2q (int nf, double as, double x, double m, Order matched_to_fixed_order) {
    check_nf(nf);
    return sxD[nf-3]->deltaMC2g(as, x, matched_to_fixed_order) * CF/CA;
  }
  //double HELLx::deltaMC2q (int nf, double as, double x, double m, Order matched_to_fixed_order) {
  //  return deltaMC2g(as, x, m, matched_to_fixed_order) * CF/CA;
  //}

  double HELLx::deltaMCLg (int nf, double as, double x, double m, Order matched_to_fixed_order) {
    check_nf(nf);
    return sxD[nf-3]->deltaMCLg(as, x, matched_to_fixed_order);
  }
  double HELLx::deltaMCLq (int nf, double as, double x, double m, Order matched_to_fixed_order) {
    check_nf(nf);
    return sxD[nf-3]->deltaMCLg(as, x, matched_to_fixed_order) * CF/CA;
  }
  //double HELLx::deltaMCLq (int nf, double as, double x, double m, Order matched_to_fixed_order) {
  //  return deltaMCLg(as, x, m, matched_to_fixed_order) * CF/CA;
  //}












};

