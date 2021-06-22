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
#include <utility>
#include "HELL/include/math/matrix.hh"
#include "HELL/include/version.hh"


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

  double Qofalphas(double as, int nf);


  // mother class for xTable files
  class xTable {
  protected:
    bool quiet;
    ifstream *infile;
    string vers;
    double *xx, *logx;
    int Np1, Np2;
    double x_min, x_mid, x_max;
    double interpolate(double x);
  public:
    xTable(string filename);
    ~xTable() { delete[] xx; delete[] logx; /*delete infile;*/ };
    void SetQuietMode(bool qu=true) { quiet = qu; }
    string GetVersion() { return vers; }
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

  // daugther class for massless DIS coefficient functions
  class xTableC : public xTable {
  private:
    double *xdC2g, *xdCLg;
    void Init();
  public:
    xTableC(string filename, bool) : xTable(filename) { Init(); }
    ~xTableC() { delete[] xdC2g; delete[] xdCLg; };
    void eval(double x, double &dC2g, double &dCLg);
  };

  // daugther class for massive DIS coefficient functions
  class xTableCm : public xTable {
  private:
    double *mQvals, **xdKhg, **xdC2g, **xdCLg, **xdC2axg, **xdCLaxg, **xdC2CCg, **xdCLCCg, **xdC3CCg;
    int Nmass;
    double Q;
    void Init();
  public:
    xTableCm(string filename, bool) : xTable(filename) { Init(); }
    ~xTableCm() { for(int i=0; i<Nmass; i++) {delete[] xdCLg[i]; delete[] xdC2g[i]; delete[] xdCLaxg[i]; delete[] xdC2axg[i]; delete[] xdCLCCg[i]; delete[] xdC2CCg[i]; delete[] xdC3CCg[i]; delete[] xdKhg[i];} delete[] xdCLg; delete[] xdC2g; delete[] xdKhg; delete[] mQvals; };
    void eval(double x, double mQ, double &dKhg, double &dC2g, double &dCLg, double &dC2axg, double &dCLaxg, double &dC2CCg, double &dCLCCg, double &dC3CCg, double as, int nf);
  };

  // daugther class for Higgs coefficient function
  class xTableCggH : public xTable {
  private:
    double *xdCggH, *xdCggHaux;
    double c10, c20, c11, c30, c21;
    void Init();
  public:
    xTableCggH(string filename, bool) : xTable(filename) { Init(); }
    ~xTableCggH() { delete[] xdCggH; delete[] xdCggHaux; };
    void eval(double x, double &dCggH, double &dCggHaux);
    double GetCoeff10() { return c10; };
    double GetCoeff20() { return c20; };
    double GetCoeff11() { return c11; };
    double GetCoeff30() { return c30; };
    double GetCoeff21() { return c21; };
  };




  // Class for nf fixed.
  //
  class HELLxnf {
  private:
    int _nf;
    int _order;
    int _RCvar;
    bool useLLp;
    bool quiet;
    vector<double> _alphas, _alphasHgg;
    string datapath_;
    map<int,xTableP*>    xT[3];
    map<int,xTableC*>    xTC[3];
    map<int,xTableCm*>   xTCm[3];
    map<int,xTableCggH*> xTCggH[3];
    int alphas_interpolation(double as, vector<double> vas, double &factor);
    double alphas_cubicinterpolate(double as, double k, vector<double> vas, double *y);
    template<class S> void ReadTable (int k, map<int,S*> &T, string basename);
    double DeltaC   (double as, double x, Order matched_to_fixed_order, double muFrat, string id);
    double DeltaCm  (double as, double x, Order matched_to_fixed_order, double muFrat, string id, double mQ);
    double DeltaCggH(double as, double x, Order matched_to_fixed_order, double muFrat, bool aux=false);
  public:
    HELLxnf(int nf, LogOrder order, string datapath="./data/") : _nf(nf), _order(order), _RCvar(0), useLLp(false), quiet(false) { Init(datapath); }
    ~HELLxnf();
    //
    void Init(string datapath);
    //
    // RCvar=0   central
    // RCvar=1   variation of gamma_plus (using LO alpha_s approximation for RC, "Airy")
    // RCvar=2   variation of gamma_qg   (using r -> as*b0)
    void SetRCvar(int var=1) { _RCvar = var; }
    int  GetOrder();
    int  GetNf();
    void GetAvailableAlphas(vector<double> &as);
    void SetLLpMode(bool use_LLp=true) { useLLp = use_LLp; }
    void SetQuietMode(bool qu=true) { quiet = qu; }
    //
    // Delta P   (LL can be matched to LO or NLO)  (NLL can be matched to NLO or NNLO)
    sqmatrix<double> DeltaP(double as, double x, Order matched_to_fixed_order = NLO);
    //
    // Delta matching condition K_hi (i=g,q)   (can be matched to NLO or NNLO)
    double deltaKhg  (double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO);
    double deltaKhq  (double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO);
    //
    // Delta coefficient functions   (can be matched to NLO or NNLO)
    // Massless DIS
    double deltaC2g  (double as, double x, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS F2
    double deltaC2q  (double as, double x, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS F2
    double deltaCLg  (double as, double x, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS FL
    double deltaCLq  (double as, double x, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS FL
    // Massive DIS, Neutral Current
    double deltaMC2g  (double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive F2
    double deltaMC2q  (double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive F2
    double deltaMCLg  (double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive FL
    double deltaMCLq  (double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive FL
    double deltaMC2axg(double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive F2, axial contribution
    double deltaMC2axq(double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive F2, axial contribution
    double deltaMCLaxg(double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive FL, axial contribution
    double deltaMCLaxq(double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive FL, axial contribution
    // Please note that massless coefficient functions contain a factor of nf in them.
    // FONLL=SACOT is obtained (for a single quark flavour, e.g. charm) as
    //   deltaC2_FONLL = deltaMC2i - deltaKhi      (i=g,q) for F2
    //   deltaCL_FONLL = deltaMCLi                 (i=g,q) for FL
    // In the m/Q->0 limit, these expressions tend to deltaCki/nf  (k=2,L i=g,q)
    //
    // Massive DIS, Charged Current
    double deltaMC2CCg(double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive F2
    double deltaMC2CCq(double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive F2
    double deltaMCLCCg(double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive FL
    double deltaMCLCCq(double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive FL
    double deltaMC3CCg(double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive F3  (this is for production of heavy quark. Heavy antiquark have an overall opposite sign)
    double deltaMC3CCq(double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive F3  (this is for production of heavy quark. Heavy antiquark have an overall opposite sign)
    // FONLL=SACOT is obtained (for a single quark flavour, e.g. charm) as
    //   deltaC2CC_FONLL = deltaMC2i - deltaKhi/2      (i=g,q) for F2
    //   deltaCLCC_FONLL = deltaMCLi                   (i=g,q) for FL
    //   deltaC2CC_FONLL = deltaMC2i + deltaKhi/2      (i=g,q) for F3  (here the sign is for production of heavy quark. Heavy antiquark have an overall opposite sign)
    // In the m/Q->0 limit, these expressions tend to deltaCki/nf  (k=2,L,3 i=g,q); note that massles C3 is zero


    double deltaCggH   (double as, double x, Order matched_to_fixed_order = NLO, double muFrat=1);  // ggH, gg channel
    double deltaCggHaux(double as, double x, Order matched_to_fixed_order = NLO, double muFrat=1);  // ggH, auxiliary integral for qg and qq channels
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
    sqmatrix<double>   DeltaP(int nf, double as, double x, Order matched_to_fixed_order = NLO);
    //
    // Delta matching condition K_hi (i=g,q)
    double deltaKhg  (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO);
    double deltaKhq  (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO);
    //
    // Delta coefficient functions
    // Massless DIS
    double deltaC2g  (int nf, double as, double x, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS F2
    double deltaC2q  (int nf, double as, double x, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS F2
    double deltaCLg  (int nf, double as, double x, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS FL
    double deltaCLq  (int nf, double as, double x, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS FL
    // Massive DIS, Neutral Current
    double deltaMC2g  (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive F2
    double deltaMC2q  (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive F2
    double deltaMCLg  (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive FL
    double deltaMCLq  (int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive FL
    double deltaMC2axg(int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive F2, axial contribution
    double deltaMC2axq(int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive F2, axial contribution
    double deltaMCLaxg(int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive FL, axial contribution
    double deltaMCLaxq(int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive FL, axial contribution
    // Please note that massless coefficient functions contain a factor of nf in them.
    // FONLL=SACOT is obtained (for a single quark flavour, e.g. charm) as
    //   deltaC2_FONLL = deltaMC2i - deltaKhi      (i=g,q) for F2
    //   deltaCL_FONLL = deltaMCLi                 (i=g,q) for FL
    // In the m/Q->0 limit, these expressions tend to deltaCki/nf  (k=2,L i=g,q)
    //
    // Massive DIS, Charged Current
    double deltaMC2CCg(int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive F2
    double deltaMC2CCq(int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive F2
    double deltaMCLCCg(int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive FL
    double deltaMCLCCq(int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive FL
    double deltaMC3CCg(int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive F3  (this is for production of heavy quark. Heavy antiquark have an overall opposite sign)
    double deltaMC3CCq(int nf, double as, double x, double m_Q_ratio, Order matched_to_fixed_order = NLO, double muFrat=1);  // DIS massive F3  (this is for production of heavy quark. Heavy antiquark have an overall opposite sign)
    // FONLL=SACOT is obtained (for a single quark flavour, e.g. charm) as
    //   deltaC2CC_FONLL = deltaMC2i - deltaKhi/2      (i=g,q) for F2
    //   deltaCLCC_FONLL = deltaMCLi                   (i=g,q) for FL
    //   deltaC2CC_FONLL = deltaMC2i + deltaKhi/2      (i=g,q) for F3  (here the sign is for production of heavy quark. Heavy antiquark have an overall opposite sign)
    // In the m/Q->0 limit, these expressions tend to deltaCki/nf  (k=2,L,3 i=g,q); note that massles C3 is zero
    //
  };




};



