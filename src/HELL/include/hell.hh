/* -----------------------------------------

   _  _ ____ _    _    
   |__| |___ |    |    
   |  | |___ |___ |___ x
   
   HELLx: High-Energy Large Logarithms - fast x-space version

   Author: Marco Bonvini

   Computes Delta P_res = P_res - P_FixedOrder
   as an interpolation from pre-prepared input files

   ----------------------------------------- */

#pragma once

#include <iostream>
#include <vector>
#include <map>
#include "./math/matrix.hh"


using namespace std;


#define gg entry11
#define gq entry12
#define qg entry21
#define qq entry22



namespace HELLx {

  const double CA = 3.;
  const double CF = 0.5*CA - 0.5/CA;
  enum Order {
    none = -1,
    LO = 0,
    NLO = 1,
    NNLO = 2,
    N2LO = 2,
    N3LO  = 3,
    NNNLO = 3
  };
  enum LogOrder {
    LL = 0,
    NLL = 1
  };


  // mother class for xTable files
  class xTable {
  protected:
    ifstream *infile;
    double *xx;
    int Np1, Np2;
    double x_min, x_mid, x_max;
    double interpolate(double x);
  public:
    xTable(string filename);
    //~xTable() { delete[] xx; delete infile; };
    ~xTable() { delete[] xx; };
  };

  // daugther class for splitting functions
  class xTableP : public xTable {
  private:
    bool isNLL;
    double *xdPgg, *xdPqg;
    void Init();
  public:
    xTableP(string filename, bool nll) : xTable(filename) { isNLL=nll; Init(); }
    ~xTableP() { delete[] xdPgg; delete[] xdPqg; };
    void eval(double x, double &dPgg, double &dPqg);
  };

  // daugther class for massless coefficient functions
  class xTableC : public xTable {
  private:
    double *xdC2g, *xdCLg;
    void Init();
  public:
    xTableC(string filename, bool nll) : xTable(filename) { Init(); }
    ~xTableC() { delete[] xdC2g; delete[] xdCLg; };
    void eval(double x, double &dC2g, double &dCLg);
  };

  // daugther class for massive coefficient functions
  class xTableCm : public xTable {
  private:
    double *mQvals, **xdKhg, **xdC2g, **xdCLg;
    int Nmass;
    double Q;
    void Init();
  public:
    xTableCm(string filename, bool nll) : xTable(filename) { Init(); }
    ~xTableCm() { for(int i=0; i<Nmass; i++) {delete[] xdCLg[i]; delete[] xdC2g[i]; delete[] xdKhg[i];} delete[] xdCLg; delete[] xdC2g; delete[] xdKhg; delete[] mQvals; };
    void eval(double x, double mQ, double &dKhg, double &dC2g, double &dCLg, double as, int nf);
  };




  // Class for nf fixed.
  //
  class HELLxnf {
  private:
    int _nf;
    int _order;
    vector<double> _alphas;
    string datapath_;
    map<int,xTableP*>  xT;
    map<int,xTableC*>  xTC;
    map<int,xTableCm*> xTCm;
    //map<int,xTable*>::iterator itxT;
    int alphas_interpolation(double as, vector<double> vas, double &factor);
    template<class S> void ReadTable (int k, map<int,S*> &T, string basename);
    double DeltaC (double as, double x, Order matched_to_fixed_order, string id);
    double DeltaCm(double as, double x, Order matched_to_fixed_order, string id, double mQ);
  public:
    HELLxnf(int nf, LogOrder order, string datapath="./data/") : _nf(nf), _order(order) { Init(datapath); }
    ~HELLxnf();
    //
    void Init(string datapath);
    //
    int GetOrder();
    int GetNf();
    void GetAvailableAlphas(vector<double> &as);
    //
    // Delta P   (LL can be matched to LO or NLO)  (NLL can be matched to NLO or NNLO)
    sqmatrix<double> DeltaP(double as, double x,  Order matched_to_fixed_order = NLO);
    //
    // Delta matching condition K_hi (i=g,q)   (can be matched to NLO or NNLO)
    double deltaKhg  (double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO);
    double deltaKhq  (double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO);
    //
    // Delta coefficient functions   (can be matched to NLO or NNLO)
    double deltaC2g  (double as, double x, Order matched_to_fixed_order = NLO);  // DIS F2
    double deltaC2q  (double as, double x, Order matched_to_fixed_order = NLO);  // DIS F2
    double deltaCLg  (double as, double x, Order matched_to_fixed_order = NLO);  // DIS FL
    double deltaCLq  (double as, double x, Order matched_to_fixed_order = NLO);  // DIS FL
    double deltaMC2g (double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO);  // DIS massive F2
    double deltaMC2q (double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO);  // DIS massive F2
    double deltaMCLg (double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO);  // DIS massive FL
    double deltaMCLq (double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO);  // DIS massive FL
    // Please note that massless coefficient functions contain a factor of nf in them.
    // FONLL=SACOT is obtained (for a single quark flavour, e.g. charm) as
    //   deltaC2_FONLL = deltaMC2i - deltaKhi      (i=g,q) for F2
    //   deltaCL_FONLL = deltaMCLi                 (i=g,q) for FL
    // In the m/Q->0 limit, these expressions tend to deltaCki/nf  (k=2,L i=g,q)
  };





  // Class for variable nf dependence
  //
  class HELLx {
  private:
    HELLxnf *sxD[4]; // one for each nf
  public:
    HELLx(LogOrder order, string datapath="./data/");
    ~HELLx();
    //
    HELLxnf* GetHELLxnf(int nf);
    //
    // Delta P
    sqmatrix<double>   DeltaP(int nf, double as, double x,  Order matched_to_fixed_order = NLO);
    //
    // Delta matching condition K_hi (i=g,q)
    double deltaKhg  (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO);
    double deltaKhq  (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO);
    //
    // Delta coefficient functions
    double deltaC2g  (int nf, double as, double x, Order matched_to_fixed_order = NLO);  // DIS F2
    double deltaC2q  (int nf, double as, double x, Order matched_to_fixed_order = NLO);  // DIS F2
    double deltaCLg  (int nf, double as, double x, Order matched_to_fixed_order = NLO);  // DIS FL
    double deltaCLq  (int nf, double as, double x, Order matched_to_fixed_order = NLO);  // DIS FL
    double deltaMC2g (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO);  // DIS massive F2
    double deltaMC2q (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO);  // DIS massive F2
    double deltaMCLg (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO);  // DIS massive FL
    double deltaMCLq (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO);  // DIS massive FL
    // Please note that massless coefficient functions contain a factor of nf in them.
    // FONLL=SACOT is obtained (for a single quark flavour, e.g. charm) as
    //   deltaC2_FONLL = deltaMC2i - deltaKhi      (i=g,q) for F2
    //   deltaCL_FONLL = deltaMCLi                 (i=g,q) for FL
    // In the m/Q->0 limit, these expressions tend to deltaCki/nf  (k=2,L i=g,q)
  };



};



