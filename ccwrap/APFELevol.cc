// C++ definition

#include "APFEL/APFELevol.h"
#include "APFEL/APFELfwevol.h"

#include <sstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
using namespace std;

#define STR_EXPAND(top) #top
#define STR(tok) STR_EXPAND(tok)

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
    std::vector<char> cstr(fname.c_str(), fname.c_str() + fname.size() + 1);
    return fexternalevolutionoperator(cstr.data(),&i,&j,&x,&beta,cstr.size()-1);
  }

  double ExternalEvolutionMatrixEv2Ev(int i, int j, int alpha, int beta)
  {
    return fexternalevolutionmatrixev2ev(&i,&j,&alpha,&beta);
  }

  double ExternalEvolutionMatrixEv2Ph(int i, int j, int alpha, int beta)
  {
    return fexternalevolutionmatrixev2ph(&i,&j,&alpha,&beta);
  }

  double ExternalEvolutionMatrixPh2Ph(int i, int j, int alpha, int beta)
  {
    return fexternalevolutionmatrixph2ph(&i,&j,&alpha,&beta);
  }

  void ComputeExternalSplittingFunctions(const std::string& fname, int pt, int nf, double x, int beta)
  {
    std::vector<char> cstr(fname.c_str(), fname.c_str() + fname.size() + 1);
    return fcomputeexternalsplittingfunctions(cstr.data(),&pt,&nf,&x,&beta,cstr.size()-1);
  }

  double ExternalSplittingFunctions(int i, int j)
  {
    return fexternalsplittingfunctions(&i,&j);
  }

  void LHAPDFgrid(int Nrep, double Qin, const std::string& fname)
  {
    std::vector<char> cstr(fname.c_str(), fname.c_str() + fname.size() + 1);
    flhapdfgrid(&Nrep, &Qin, cstr.data(), cstr.size()-1);
  }

  void LHAPDFgridDerivative(int Nrep, const std::string& fname)
  {
    std::vector<char> cstr(fname.c_str(), fname.c_str() + fname.size() + 1);
    flhapdfgridderivative(&Nrep, cstr.data(), cstr.size()-1);
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
    std::vector<char> cstr(la.c_str(), la.c_str() + la.size() + 1);
    fsetsmallxresummation(&sx, cstr.data(), cstr.size()-1);
  }

  double HeavyQuarkMass(int i,double Q)
  {
    return fheavyquarkmass(&i,&Q);
  }

  double GetThreshold(int i)
  {
    return fgetthreshold(&i);
  }

  int GetMaxFlavourAlpha()
  {
    return fgetmaxflavouralpha();
  }

  int GetMaxFlavourPDFs()
  {
    return fgetmaxflavourpdfs();
  }

  double HeavyQuarkThreshold(int i)
  {
    return fheavyquarkthreshold(&i);
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
    std::vector<char> cstr(evol.c_str(), evol.c_str() + evol.size() + 1);
    fsetalphaevolution(cstr.data(), cstr.size()-1);
  }

  void SetLambdaQCDRef(double lambdaref, int nref)
  {
    fsetlambdaqcdref(&lambdaref,&nref);
  }

  void SetPDFEvolution(const std::string& evolp)
  {
    std::vector<char> cstr(evolp.c_str(), evolp.c_str() + evolp.size() + 1);
    fsetpdfevolution(cstr.data(), cstr.size()-1);
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
    std::vector<char> cstr(name.c_str(), name.c_str() + name.size() + 1);
    fsetpdfset(cstr.data(), cstr.size()-1);
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
    std::vector<char> cstr(theory.c_str(), theory.c_str() + theory.size() + 1);
    fsettheory(cstr.data(), cstr.size()-1);
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
