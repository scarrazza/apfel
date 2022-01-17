// C++ definition

#include "APFEL/APFELobs.h"
#include "APFEL/APFELfwobs.h"

#include <sstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
using namespace std;

#define STR_EXPAND(top) #top
#define STR(tok) STR_EXPAND(tok)

namespace APFEL {

  void InitializeAPFEL_DIS()
  {
    finitializeapfel_dis();
  }

  void ComputeStructureFunctionsAPFEL(double Q0, double Q)
  {
    fcomputestructurefunctionsapfel(&Q0,&Q);
  }

  void CacheStructureFunctionsAPFEL(double Q0)
  {
    fcachestructurefunctionsapfel(&Q0);
  }

  void SetMassScheme(const std::string& ms)
  {
    std::vector<char> cms(ms.c_str(), ms.c_str() + ms.size() + 1);
    fsetmassscheme(cms.data());
  }

  void SetPolarizationDIS(double pol)
  {
    fsetpolarizationdis(&pol);
  }

  void SetProcessDIS(const std::string& pr)
  {
    std::vector<char> cpr(pr.c_str(), pr.c_str() + pr.size() + 1);
    fsetprocessdis(cpr.data());
  }

  void SetProjectileDIS(const std::string& lept)
  {
    std::vector<char> clept(lept.c_str(), lept.c_str() + lept.size() + 1);
    fsetprojectiledis(clept.data());
  }

  void SetTargetDIS(const std::string& tar)
  {
    std::vector<char> ctar(tar.c_str(), tar.c_str() + tar.size() + 1);
    fsettargetdis(ctar.data());
  }

  void SelectCharge(const std::string& selch)
  {
    std::vector<char> cselch(selch.c_str(), selch.c_str() + selch.size() + 1);
    fselectcharge(cselch.data());
  }

  double ExternalDISOperator(const std::string& SF,int ihq,int i,double x,int beta)
  {
    std::vector<char> cSF(SF.c_str(), SF.c_str() + SF.size() + 1);
    return fexternaldisoperator(cSF.data(),&ihq,&i,&x,&beta);
  }

  double F2light(double x)
  {
    return ff2light(&x);
  }

  double F2charm(double x)
  {
    return ff2charm(&x);
  }

  double F2bottom(double x)
  {
    return ff2bottom(&x);
  }

  double F2top(double x)
  {
    return ff2top(&x);
  }

  double F2total(double x)
  {
    return ff2total(&x);
  }

  double FLlight(double x)
  {
    return ffllight(&x);
  }

  double FLcharm(double x)
  {
    return fflcharm(&x);
  }

  double FLbottom(double x)
  {
    return fflbottom(&x);
  }

  double FLtop(double x)
  {
    return ffltop(&x);
  }

  double FLtotal(double x)
  {
    return ffltotal(&x);
  }

  double F3light(double x)
  {
    return ff3light(&x);
  }

  double F3charm(double x)
  {
    return ff3charm(&x);
  }

  double F3bottom(double x)
  {
    return ff3bottom(&x);
  }

  double F3top(double x)
  {
    return ff3top(&x);
  }

  double F3total(double x)
  {
    return ff3total(&x);
  }

  double g1light(double x)
  {
    return fg1light(&x);
  }

  double g1charm(double x)
  {
    return fg1charm(&x);
  }

  double g1bottom(double x)
  {
    return fg1bottom(&x);
  }

  double g1top(double x)
  {
    return fg1top(&x);
  }

  double g1total(double x)
  {
    return fg1total(&x);
  }

  double gLlight(double x)
  {
    return fgllight(&x);
  }

  double gLcharm(double x)
  {
    return fglcharm(&x);
  }

  double gLbottom(double x)
  {
    return fglbottom(&x);
  }

  double gLtop(double x)
  {
    return fgltop(&x);
  }

  double gLtotal(double x)
  {
    return fgltotal(&x);
  }

  double g4light(double x)
  {
    return fg4light(&x);
  }

  double g4charm(double x)
  {
    return fg4charm(&x);
  }

  double g4bottom(double x)
  {
    return fg4bottom(&x);
  }

  double g4top(double x)
  {
    return fg4top(&x);
  }

  double g4total(double x)
  {
    return fg4total(&x);
  }

  double StructureFunctionxQ(const std::string& proc, const std::string& sf, const std::string& comp, double x, double Q)
  {
    std::vector<char> cproc(proc.c_str(), proc.c_str() + proc.size() + 1);
    std::vector<char> csf(sf.c_str(), sf.c_str() + sf.size() + 1);
    std::vector<char> ccomp(comp.c_str(), comp.c_str() + comp.size() + 1);
    return fstructurefunctionxq(cproc.data(),csf.data(),ccomp.data(),&x,&Q);
  }

  void SetZMass(double massz)
  {
    fsetzmass(&massz);
  }

  void SetWMass(double massw)
  {
    fsetwmass(&massw);
  }

  void SetProtonMass(double massp)
  {
    fsetprotonmass(&massp);
  }

  void SetSin2ThetaW(double sw)
  {
    fsetsin2thetaw(&sw);
  }

  void SetCKM(double vud,double vus,double vub,double vcd,double vcs,double vcb,double vtd,double vts,double vtb)
  {
    fsetckm(&vud,&vus,&vub,&vcd,&vcs,&vcb,&vtd,&vts,&vtb);
  }

  void SetPropagatorCorrection(double dr)
  {
    fsetpropagatorcorrection(&dr);
  }

  void SetEWCouplings(double vd,double vu,double ad,double au)
  {
    fsetewcouplings(&vd,&vu,&ad,&au);
  }

  void SetGFermi(double gf)
  {
    fsetgfermi(&gf);
  }

  void SetRenQRatio(double ratioR)
  {
    fsetrenqratio(&ratioR);
  }

  void SetFacQRatio(double ratioF)
  {
    fsetfacqratio(&ratioF);
  }

  void EnableDynamicalScaleVariations(int dsv)
  {
    fenabledynamicalscalevariations(&dsv);
  }

  void EnableIntrinsicCharm(int ic)
  {
    fenableintrinsiccharm(&ic);
  }

  double GetZMass()
  {
    return fgetzmass();
  }

  double GetWMass()
  {
    return fgetwmass();
  }

  double GetProtonMass()
  {
    return fgetprotonmass();
  }

  double GetSin2ThetaW()
  {
    return fgetsin2thetaw();
  }

  double GetCKM(int u, int d)
  {
    return fgetckm(&u,&d);
  }

  double GetGFermi()
  {
    return fgetgfermi();
  }

  double GetSIATotalCrossSection(int pto, double q, const std::string& comp)
  {
    std::vector<char> ccomp(comp.c_str(), comp.c_str() + comp.size() + 1);
    return fgetsiatotalcrosssection(&pto,&q,ccomp.data());
  }

  void EnableTargetMassCorrections(int tc)
  {
    fenabletargetmasscorrections(&tc);
  }

  void EnableDampingFONLL(int df)
  {
    fenabledampingfonll(&df);
  }

  void SetDampingPowerFONLL(int dpc, int dpb, int dpt)
  {
    fsetdampingpowerfonll(&dpc, &dpb, &dpt);
  }

  void ComputeChargesDIS(double q2, double *bq, double *dq, double *bqt)
  {
    fcomputechargesdis(&q2,bq,dq,bqt);
  }

  double F2LO(double x, double q)
  {
    double Q2 = q * q;
    double bq[7];
    double dq[7];
    double bqt[7];
    ComputeChargesDIS(Q2,bq,dq,bqt);
    double xf[13];
    xPDFxQall(x, q, xf);
    double F2 = 0;
    for (int j = 1; j <=6; j++)
      F2 += bq[j] * ( xf[6+j] + xf[6-j] );
    return F2;
  }

  double FKSimulator(double x,double q,double y,int i,int beta)
  {
    return ffksimulator(&x,&q,&y,&i,&beta);
  }

  void SetFKObservable(const std::string& obs)
  {
    std::vector<char> cobs(obs.c_str(), obs.c_str() + obs.size() + 1);
    fsetfkobservable(cobs.data());
  }

  void GetFKObservable()
  {
    fgetfkobservable();
  }

  double FKObservables(double x,double q,double y)
  {
    return ffkobservables(&x,&q,&y);
  }

  void ComputeFKTables(const std::string& inputfile, const std::string& outputpath,
		       double Q0, int* flmap)
  {
    std::vector<char> cinputfile(inputfile.c_str(), inputfile.c_str() + inputfile.size() + 1);
    std::vector<char> coutputpath(outputpath.c_str(), outputpath.c_str() + outputpath.size() + 1);
    fcomputefktables(cinputfile.data(),coutputpath.data(),&Q0,flmap);
  }

  void ComputeHardCrossSectionsDY(const std::string& inputfile, 
				  const std::string& outputfile)
  {
    std::vector<char> cinputfile(inputfile.c_str(), inputfile.c_str() + inputfile.size() + 1);
    std::vector<char> coutputfile(outputfile.c_str(), outputfile.c_str() + outputfile.size() + 1);
    fcomputehardcrosssectionsdy(cinputfile.data(),coutputfile.data());
  }

  void EnableSFNLOQEDCorrections(int qedsfc)
  {
    fenablesfnloqedcorrections(&qedsfc);
  }

  void LHAPDFgridStructureFunctions(int Nrep, double Qin, const std::string& fname)
  {
    std::vector<char> cfname(fname.c_str(), fname.c_str() + fname.size() + 1);
    flhapdfgridstructurefunctions(&Nrep,&Qin,cfname.data());
  }

  void SetScaleVariationProcedure(int svp)
  {
    fsetscalevariationprocedure(&svp);
  }
}
