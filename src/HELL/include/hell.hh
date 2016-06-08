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
#include "HELL/include/math/matrix.hh"


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
    NNLO = 2
  };
  enum LogOrder {
    LL = 0,
    NLL = 1
  };



  class xTable {
  private:
    double *xx;
    double *xdPgg, *xdPqg;//, *xdC2g, *xdCLg;//, *xdMC2g, *xdMCLg;
    int Np1, Np2;
    double x_min, x_mid, x_max;
    bool isNLL;
  public:
    xTable(string filename, bool nll);
    ~xTable() {};
    double eval (double x);
    void   evalP(double x, double &dPgg, double &dPqg);
    //void   evalC(double x, double &dC2g, double &dCLg);
  };



  // Class for nf fixed.
  //
  class HELLxnf {
  private:
    int _nf;
    int _order;
    vector<double> _alphas;
    string datapath_;
    map<int,xTable*> xT, xTC;
    map<int,xTable*>::iterator itxT;
    int alphas_interpolation(double as, vector<double> vas, double &factor);
    void ReadTable (int k, map<int,xTable*> &T, string basename);
    double DeltaC(double as, double x, Order matched_to_fixed_order, string id);
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
    // Delta P
    sqmatrix<double> DeltaP(double as, double x,  Order matched_to_fixed_order = NLO);
    //
    // Delta Coefficient functions (the quark coefficients are found multiplying by CF/CA)
    double deltaC2g  (double as, double x, Order matched_to_fixed_order = NLO);  // DIS F2
    double deltaC2q  (double as, double x, Order matched_to_fixed_order = NLO);  // DIS F2
    double deltaCLg  (double as, double x, Order matched_to_fixed_order = NLO);  // DIS FL
    double deltaCLq  (double as, double x, Order matched_to_fixed_order = NLO);  // DIS FL
    double deltaMC2g (double as, double x, double m, Order matched_to_fixed_order = NLO);  // DIS massive F2
    double deltaMC2q (double as, double x, double m, Order matched_to_fixed_order = NLO);  // DIS massive F2
    double deltaMCLg (double as, double x, double m, Order matched_to_fixed_order = NLO);  // DIS massive FL
    double deltaMCLq (double as, double x, double m, Order matched_to_fixed_order = NLO);  // DIS massive FL
  };





  // Class for variable nf dependence
  //
  const int nf_min = 3;
  const int nf_max = 6;
  const int values_of_nf = nf_max - nf_min + 1;
  //
  class HELLx {
  private:
    int _order;
    double mass_c, mass_b, mass_t;
    double as_c, as_b, as_t;
    //bool as_thr_init;
    //bool check_as_thr_init();
    HELLxnf *sxD[values_of_nf];
    double as_LO (double mu);
    double as_NLO(double mu);
    double as_thr(double mu);
  public:
    HELLx(LogOrder order, string datapath="./data/");
    ~HELLx();
    //
    int GetOrder();
    //
    void init_mass_thresholds();
    void init_mass_thresholds(double mc, double mb, double mt);
    void init_as_thresholds();
    void init_as_thresholds(double as_of_mu(double));
    void init_as_thresholds(double asc, double asb, double ast);
    int nf_of_as(double as);
    //
    // Delta P
    sqmatrix<double>   DeltaP    (double as, double x,  Order matched_to_fixed_order = NLO);
    //
    // Delta Coefficient functions (the quark coefficients are found multiplying by CF/CA)
    double deltaC2g  (double as, double x, Order matched_to_fixed_order = NLO);  // DIS F2
    double deltaC2q  (double as, double x, Order matched_to_fixed_order = NLO);  // DIS F2
    double deltaCLg  (double as, double x, Order matched_to_fixed_order = NLO);  // DIS FL
    double deltaCLq  (double as, double x, Order matched_to_fixed_order = NLO);  // DIS FL
    double deltaMC2g (double as, double x, double m, Order matched_to_fixed_order = NLO);  // DIS massive F2
    double deltaMC2q (double as, double x, double m, Order matched_to_fixed_order = NLO);  // DIS massive F2
    double deltaMCLg (double as, double x, double m, Order matched_to_fixed_order = NLO);  // DIS massive FL
    double deltaMCLq (double as, double x, double m, Order matched_to_fixed_order = NLO);  // DIS massive FL
  };




  /*
  // Delta P_+
  extern double  xdeltaPplus   (double as, double  x);
  extern dcomplex deltaGammaplus(double as, dcomplex N);
  // Delta P_qg
  extern double xdeltaPqg    (double as, double  x);
  extern dcomplex deltaGammaqg(double as, dcomplex N);
  // Delta C_2g
  extern double xdeltaC2g  (double as, double  x);
  extern dcomplex deltaC2g_N(double as, dcomplex N);
  // Delta C_2g  ---  as-derivative
  extern double xdeltaC2g_deriv  (double as, double  x);
  extern dcomplex deltaC2g_deriv_N(double as, dcomplex N);

  //extern sqmatrix DeltaGammaLL (dcomplex N, double as, int nf);
  //extern sqmatrix DeltaGammaNLL(dcomplex N, double as, int nf);
  extern sqmatrix DeltaGamma(dcomplex N, double as, int nf);
  */







};



