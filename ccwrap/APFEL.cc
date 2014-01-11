// C++ definition

#include "APFEL/APFEL.h"
#include "APFEL/APFELfw.h"

#include <cstring>
#include <string>
#include <iostream>
using namespace std;

#define SIZE 999
#define VERSION "2.0.1"

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

  double xgamma(double x)
  {
    return fxgamma(&x);
  }

  void LHAPDFgrid(int Nrep, double Qin, const std::string& fname)
  {
    char cfname[SIZE+1];
    strncpy(cfname, fname.c_str(), SIZE);
    flhapdfgrid(&Nrep,&Qin,cfname);
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

  std::string GetVersion(void) 
  { 
    return VERSION; 
  }

  void CleanUp(void)
  {
    fcleanup();
  }

  void EnableWelcomeMessage(bool wc)
  {
    fenablewelcomemessage(&wc);
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
}
