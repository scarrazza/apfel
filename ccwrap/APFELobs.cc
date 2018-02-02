// C++ definition

#include "APFEL/APFELobs.h"
#include "APFEL/APFELfwobs.h"

#include <sstream>
#include <cstring>
#include <string>
#include <iostream>
using namespace std;

#define STR_EXPAND(top) #top
#define STR(tok) STR_EXPAND(tok)

#define SIZE 999

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
    char cms[SIZE+1];
    strncpy(cms, ms.c_str(), SIZE);
    fsetmassscheme(cms);
  }

  void SetPolarizationDIS(double pol)
  {
    fsetpolarizationdis(&pol);
  }

  void SetProcessDIS(const std::string& pr)
  {
    char cpr[SIZE+1];
    strncpy(cpr, pr.c_str(), SIZE);
    fsetprocessdis(cpr);
  }

  void SetProjectileDIS(const std::string& lept)
  {
    char clept[SIZE+1];
    strncpy(clept, lept.c_str(), SIZE);
    fsetprojectiledis(clept);
  }

  void SetTargetDIS(const std::string& tar)
  {
    char ctar[SIZE+1];
    strncpy(ctar, tar.c_str(), SIZE);
    fsettargetdis(ctar);
  }

  void SelectCharge(const std::string& selch)
  {
    char cselch[SIZE+1];
    strncpy(cselch, selch.c_str(), SIZE);
    fselectcharge(cselch);
  }

  double ExternalDISOperator(const std::string& SF,int ihq,int i,double x,int beta)
  {
    char cSF[SIZE+1];
    strncpy(cSF, SF.c_str(), SIZE);
    return fexternaldisoperator(cSF,&ihq,&i,&x,&beta);
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

  double StructureFunctionxQ(const std::string& proc, const std::string& sf, const std::string& comp, double x, double Q)
  {
    char cproc[SIZE+1];
    strncpy(cproc, proc.c_str(), SIZE);
    char csf[SIZE+1];
    strncpy(csf, sf.c_str(), SIZE);
    char ccomp[SIZE+1];
    strncpy(ccomp, comp.c_str(), SIZE);
    return fstructurefunctionxq(cproc,csf,ccomp,&x,&Q);
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
    char ccomp[SIZE+1];
    strncpy(ccomp, comp.c_str(), SIZE);
    return fgetsiatotalcrosssection(&pto,&q,ccomp);
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
    char cobs[SIZE+1];
    strncpy(cobs, obs.c_str(), SIZE);
    fsetfkobservable(cobs);
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
    char cinputfile[SIZE+1];
    strncpy(cinputfile, inputfile.c_str(), SIZE);
    char coutputpath[SIZE+1];
    strncpy(coutputpath, outputpath.c_str(), SIZE);
    fcomputefktables(cinputfile,coutputpath,&Q0,flmap);    
  }

  void ComputeHardCrossSectionsDY(const std::string& inputfile, 
				  const std::string& outputfile)
  {
    char cinputfile[SIZE+1];
    strncpy(cinputfile, inputfile.c_str(), SIZE);
    char coutputfile[SIZE+1];
    strncpy(coutputfile, outputfile.c_str(), SIZE);
    
    fcomputehardcrosssectionsdy(cinputfile,coutputfile);    
  }

  void EnableSFNLOQEDCorrections(int qedsfc)
  {
    fenablesfnloqedcorrections(&qedsfc);
  }

}
