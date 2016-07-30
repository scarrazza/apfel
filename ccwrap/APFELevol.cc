// C++ definition

#include "APFEL/APFELevol.h"
#include "APFEL/APFELfwevol.h"

#include <sstream>
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

  void DeriveAPFEL(double Q)
  {
    fderiveapfel(&Q);
  }

  void CachePDFsAPFEL(double Q0)
  {
    fcachepdfsapfel(&Q0);
  }

  double xPDF(int i, double x)
  {
    return fxpdf(&i, &x);
  }

  double xPDFxQ(int i, double x, double Q)
  {
    return fxpdfxq(&i, &x, &Q);
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

  void xPDFall(double x, double *xf)
  {
    fxpdfall(&x,xf);
  }

  void xPDFallPhoton(double x, double *xf)
  {
    fxpdfallphoton(&x,xf);
  }

  void xPDFxQall(double x, double Q, double *xf)
  {
    fxpdfxqall(&x,&Q,xf);
  }

  double xLepton(int i, double x)
  {
    return fxlepton(&i, &x);
  }

  double xLeptonj(int i, double x)
  {
    return fxleptonj(&i, &x);
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

  void EnableLeptonEvolution(int le)
  {
    fenableleptonevolution(&le);
  }

  void LockGrids(int lg)
  {
    flockgrids(&lg);
  }

  void SetTimeLikeEvolution(int tl)
  {
    fsettimelikeevolution(&tl);
  }

  void SetPolarizedEvolution(int polev)
  {
    fsetpolarizedevolution(&polev);
  }

  void SetFastEvolution(int fe)
  {
    fsetfastevolution(&fe);
  }

  void EnableMassRunning(int mr)
  {
    fenablemassrunning(&mr);
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

  double GetThreshold(int i)
  {
    return fgetthreshold(&i);
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

  void SetEpsilonTruncation(double eps)
  {
    fsetepsilontruncation(&eps);
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

  void SetQGridParameters(int npq, int degq)
  {
    fsetqgridparameters(&npq,&degq);
  }

  void SetLHgridParameters(int nx, int nxm, double xmin, double xm, double xmax,
			   int nq2, double q2min, double q2max)
  {
    fsetlhgridparameters(&nx,&nxm,&xmin,&xm,&xmax,&nq2,&q2min,&q2max);
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

  void SetMassScaleReference(double Qc, double Qb, double Qt)
  {
    fsetmassscalereference(&Qc,&Qb,&Qt);
  }

  void SetMassMatchingScales(double kmc, double kmb, double kmt)
  {
    fsetmassmatchingscales(&kmc,&kmb,&kmt);
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

  int GetPerturbativeOrder()
  {
    return fgetperturbativeorder();
  }

  double GetMuF()
  {
    return fgetmuf();
  }

  double GetMuF0()
  {
    return fgetmuf0();
  }

  void SetPoleMasses(double mc, double mb, double mt)
  {
    fsetpolemasses(&mc,&mb,&mt);
  }

  void SetTauMass(double masst)
  {
    fsettaumass(&masst);
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

  void EnableNLOQEDCorrections(int qedc)
  {
    fenablenloqedcorrections(&qedc);
  }
  
  void SetVFNS()
  {
    fsetvfns();
  }

  void ListFunctions(void) 
  {
    flistfunctions();
  }

  bool CheckAPFEL(void)
  {
    return fcheckapfel();
  }
}
