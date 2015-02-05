// C++ definition

#include "APFEL/APFEL.h"
#include "APFEL/APFELfw.h"

#include <cstring>
#include <string>
#include <iostream>
using namespace std;

#define STR_EXPAND(top) #top
#define STR(tok) STR_EXPAND(tok)

#define SIZE 999

namespace APFEL {

  void InitializeAPFEL(void) 
  {
    finitializeapfel();
  }

  void EvolveAPFEL(double Q0, double Q)
  {
    fevolveapfel(&Q0, &Q);
  }

  double xPDF(int i, double x)
  {
    return fxpdf(&i, &x);
  }

  double dxPDF(int i, double x)
  {
    return fdxpdf(&i, &x);
  }

  double xPDFj(int i, double x)
  {
    return fxpdfj(&i, &x);
  }

  double xgamma(double x)
  {
    return fxgamma(&x);
  }

  double xgammaj(double x)
  {
    return fxgammaj(&x);
  }

  double dxgamma(double x)
  {
    return fdxgamma(&x);
  }

  double ExternalEvolutionOperator(const std::string& fname, int i, int j, double x, int beta)
  {
    char cfname[SIZE+1];
    strncpy(cfname, fname.c_str(), SIZE);
    return fexternalevolutionoperator(cfname,&i,&j,&x,&beta);
  }

  void LHAPDFgrid(int Nrep, double Qin, const std::string& fname)
  {
    char cfname[SIZE+1];
    strncpy(cfname, fname.c_str(), SIZE);
    flhapdfgrid(&Nrep,&Qin,cfname);
  }

  void LHAPDFgridDerivative(int Nrep, const std::string& fname)
  {
    char cfname[SIZE+1];
    strncpy(cfname, fname.c_str(), SIZE);
    flhapdfgridderivative(&Nrep,cfname);
  }

  double AlphaQCD(double Q)
  {
    return falphaqcd(&Q);
  }

  double AlphaQED(double Q)
  {
    return falphaqed(&Q);
  }

  double NPDF(int i, int N)
  {
    return fnpdf(&i,&N);
  }
  
  double Ngamma(int N)
  {
    return fngamma(&N);
  }

  double LUMI(int i, int j, double S)
  {
    return flumi(&i,&j,&S);
  }

  double xGrid(int alpha)
  {
    return fxgrid(&alpha);
  }

  int nIntervals()
  {
    return fnintervals();
  }

  std::string GetVersion(void) 
  { 
    return STR(APFEL_VERSION); 
  }

  void CleanUp(void)
  {
    fcleanup();
  }

  void EnableWelcomeMessage(int wc)
  {
    fenablewelcomemessage(&wc);
  }

  void EnableEvolutionOperator(int eo)
  {
    fenableevolutionoperator(&eo);
  }

  void LockGrids(int lg)
  {
    flockgrids(&lg);
  }

  void SetTimeLikeEvolution(int tl)
  {
    fsettimelikeevolution(&tl);
  }

  void SetFastEvolution(int fe)
  {
    fsetfastevolution(&fe);
  }

  void SetSmallxResummation(int sx, const std::string& la)
  {
    char cla[SIZE+1];
    strncpy(cla, la.c_str(), SIZE);
    fsetsmallxresummation(&sx,cla);
  }

  double HeavyQuarkMass(int i,double Q)
  {
    return fheavyquarkmass(&i,&Q);
  }

  void SetAlphaQCDRef(double alpharef, double Qref)
  {
    fsetalphaqcdref(&alpharef,&Qref);
  }

  void SetAlphaQEDRef(double alpharef, double Qref)
  {
    fsetalphaqedref(&alpharef,&Qref);
  }

  void SetAlphaEvolution(const std::string& evol)
  {
    char cevol[SIZE+1];
    strncpy(cevol, evol.c_str(), SIZE);
    fsetalphaevolution(cevol);
  }

  void SetLambdaQCDRef(double lambdaref, int nref)
  {
    fsetlambdaqcdref(&lambdaref,&nref);
  }

  void SetPDFEvolution(const std::string& evolp)
  {
    char cevolp[SIZE+1];
    strncpy(cevolp, evolp.c_str(), SIZE);
    fsetpdfevolution(cevolp);
  }

  void SetQLimits(double Qmin, double Qmax)
  {
    fsetqlimits(&Qmin, &Qmax);
  }

  void SetFFNS(int nfl)
  {
    fsetffns(&nfl);
  }

  void SetGridParameters(int i, int np, int deg, double x)
  {
    fsetgridparameters(&i,&np,&deg,&x);
  }

  void SetExternalGrid(int i, int np, int deg, double *x)
  {
    fsetexternalgrid(&i,&np,&deg,x);
  }

  void SetMaxFlavourAlpha(int nf)
  {
    fsetmaxflavouralpha(&nf);
  }

  void SetMaxFlavourPDFs(int nf)
  {
    fsetmaxflavourpdfs(&nf);
  }

  void SetMSbarMasses(double mc, double mb, double mt)
  {
    fsetmsbarmasses(&mc,&mb,&mt);
  }
  
  void SetNumberOfGrids(int n)
  {
    fsetnumberofgrids(&n);
  }

  void SetPDFSet(const std::string& name)
  {
    char cname[SIZE+1];
    strncpy(cname, name.c_str(), SIZE);
    fsetpdfset(cname);
  }

  void SetPerturbativeOrder(int pto)
  {
    fsetperturbativeorder(&pto);
  }

  void SetPoleMasses(double mc, double mb, double mt)
  {
    fsetpolemasses(&mc,&mb,&mt);
  }

  void SetRenFacRatio(double ratio)
  {
    fsetrenfacratio(&ratio);
  }
  
  void SetReplica(int nr)
  {
    fsetreplica(&nr);
  }
  
  void SetTheory(const std::string& theory)
  {
    char ctheory[SIZE+1];
    strncpy(ctheory, theory.c_str(), SIZE);
    fsettheory(ctheory);
  }
  
  void SetVFNS()
  {
    fsetvfns();
  }
  
  void DIS_xsec(double x,double qi,double qf,double y,double pol,
		const std::string& proc,const std::string& scheme,
		int pto,const std::string& pdfset, int irep,
		const std::string& target, const std::string& proj,
		double *F2, double *F3, double *FL, double *sigma)
  {
    char cproc[SIZE+1];
    strncpy(cproc, proc.c_str(), SIZE);

    char cscheme[SIZE+1];
    strncpy(cscheme, scheme.c_str(), SIZE);

    char cpdfset[SIZE+1];
    strncpy(cpdfset, pdfset.c_str(), SIZE);

    char ctarget[SIZE+1];
    strncpy(ctarget, target.c_str(), SIZE);    

    char cproj[SIZE+1];
    strncpy(cproj, proj.c_str(), SIZE);

    fdisxsec(&x,&qi,&qf,&y,&pol,cproc,cscheme,&pto,cpdfset,&irep,ctarget, cproj,
	     F2,F3,FL,sigma);
  }

  // Functions for the new DIS module
  void InitializeAPFEL_DIS()
  {
    finitializeapfel_dis();
  }

  void ComputeStructureFunctionsAPFEL(double Q0, double Q)
  {
    fcomputestructurefunctionsapfel(&Q0,&Q);
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

  void SetDinThetaW(double sw)
  {
    fsetsinthetaw(&sw);
  }

  void SetCKM(double vud,double vus,double vub,double vcd,double vcs,double vcb,double vtd,double vts,double vtb)
  {
    fsetckm(&vud,&vus,&vub,&vcd,&vcs,&vcb,&vtd,&vts,&vtb);
  }

  void SetGFermi(double gf)
  {
    fsetgfermi(&gf);
  }

  double FKSimulator(const std::string& obs,double x,double q,double y,int i,int beta)
  {
    char cobs[SIZE+1];
    strncpy(cobs, obs.c_str(), SIZE);
    return ffksimulator(cobs,&x,&q,&y,&i,&beta);
  }

}
