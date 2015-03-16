#ifndef APFELfw_H
#define APFELfw_H

// Declarations for the Fortran/C interface

#include "APFEL/FortranWrappers.h"

extern "C" {
  
#define finitializeapfel FC_FUNC(initializeapfel,INITIALIZEAPFEL)
  void finitializeapfel(void);

#define fevolveapfel FC_FUNC(evolveapfel,EVOLVEAPFEL)
  void fevolveapfel(double*,double*);

#define fderiveapfel FC_FUNC(deriveapfel,DERIVEAPFEL)
  void fderiveapfel(double*);

#define fxpdf FC_FUNC(xpdf,XPDF)
  double fxpdf(int*, double*);

#define fxpdfj FC_FUNC(xpdfj,XPDFJ)
  double fxpdfj(int*, double*);

#define fdxpdf FC_FUNC(dxpdf,DXPDF)
  double fdxpdf(int*, double*);

#define fxgamma FC_FUNC(xgamma,XGAMMA)
  double fxgamma(double*);

#define fxgammaj FC_FUNC(xgammaj,XGAMMAJ)
  double fxgammaj(double*);

#define fdxgamma FC_FUNC(dxgamma,DXGAMMA)
  double fdxgamma(double*);

#define fxpdfall FC_FUNC(xpdfall,XPDFALL)
  void fxpdfall(double*,double*);

#define fexternalevolutionoperator FC_FUNC(externalevolutionoperator,EXTERNALEVOLUTIONOPERATOR)
  double fexternalevolutionoperator(char*,int*,int*,double*,int*);

#define flhapdfgrid FC_FUNC(lhapdfgrid,LHAPDFGRID)
  void flhapdfgrid(int*, double*, char*);

#define flhapdfgridderivative FC_FUNC(lhapdfgridderivative,LHAPDFGRIDDERIVATIVE)
  void flhapdfgridderivative(int*, char*);

#define falphaqcd FC_FUNC(alphaqcd,ALPHAQCD)
  double falphaqcd(double*);

#define falphaqed FC_FUNC(alphaqed,ALPHAQED)
  double falphaqed(double*);

#define fnpdf FC_FUNC(npdf,NPDF)
  double fnpdf(int*,int*);

#define fngamma FC_FUNC(ngamma,NGAMMA)
  double fngamma(int*);

#define flumi FC_FUNC(lumi,LUMI)
  double flumi(int*,int*,double*);

#define fxgrid FC_FUNC(xgrid,XGRID)
  double fxgrid(int*);

#define fnintervals FC_FUNC(nintervals,NINTERVALS)
  int fnintervals();

#define fcleanup FC_FUNC(cleanup,CLEANUP)
  void fcleanup(void);

#define fenablewelcomemessage FC_FUNC(enablewelcomemessage,ENABLEWELCOMEMESSAGE)
  void fenablewelcomemessage(int*);

#define fenableevolutionoperator FC_FUNC(enableevolutionoperator,ENABLEEVOLUTIONOPERATOR)
  void fenableevolutionoperator(int*);

#define flockgrids FC_FUNC(lockgrids,LOCKGRIDS)
  void flockgrids(int*);

#define fsettimelikeevolution FC_FUNC(settimelikeevolution,SETTIMELIKEEVOLUTION)
  void fsettimelikeevolution(int*);

#define fsetfastevolution FC_FUNC(setfastevolution,SETFASTEVOLUTION)
  void fsetfastevolution(int*);

#define fenablemassrunning FC_FUNC(enablemassrunning,ENABLEMASSRUNNING)
  void fenablemassrunning(int*);

#define fsetsmallxresummation FC_FUNC(setsmallxresummation,SETSMALLXRESUMMATION)
  void fsetsmallxresummation(int*,char*);

#define fheavyquarkmass FC_FUNC(heavyquarkmass,HEAVYQUARKMASS)
  double fheavyquarkmass(int*,double*);

#define fsetalphaqcdref FC_FUNC(setalphaqcdref,SETALPHAQCDREF)
  void fsetalphaqcdref(double*,double*);

#define fsetalphaqedref FC_FUNC(setalphaqedref,SETALPHAQEDREF)
  void fsetalphaqedref(double*,double*);

#define fsetalphaevolution FC_FUNC(setalphaevolution,SETALPHAEVOLUTION)
  void fsetalphaevolution(char*);

#define fsetlambdaqcdref FC_FUNC(setlambdaqcdref,SETLAMBDAQCDREF)
  void fsetlambdaqcdref(double*,int*);

#define fsetpdfevolution FC_FUNC(setpdfevolution,SETPDFEVOLUTION)
  void fsetpdfevolution(char*);

#define fsetqlimits FC_FUNC(setqlimits,SETQLIMITS)
  void fsetqlimits(double*,double*);

#define fsetffns FC_FUNC(setffns,SETFFNS)
  void fsetffns(int*);

#define fsetgridparameters FC_FUNC(setgridparameters,SETGRIDPARAMETERS)
  void fsetgridparameters(int*,int*,int*,double*);

#define fsetexternalgrid FC_FUNC(setexternalgrid,SETEXTERNALGRID)
  void fsetexternalgrid(int*,int*,int*,double*);

#define fsetmaxflavouralpha FC_FUNC(setmaxflavouralpha,SETMAXFLAVOURALPHA)
  void fsetmaxflavouralpha(int*);

#define fsetmaxflavourpdfs FC_FUNC(setmaxflavourpdfs,SETMAXFLAVOURPDFS)
  void fsetmaxflavourpdfs(int*);

#define fsetmsbarmasses FC_FUNC(setmsbarmasses,SETMSBARMASSES)
  void fsetmsbarmasses(double*,double*,double*);

#define fsetnumberofgrids FC_FUNC(setnumberofgrids,SETNUMBEROFGRIDS)
  void fsetnumberofgrids(int*);

#define fsetpdfset FC_FUNC(setpdfset,SETPDFSET)
  void fsetpdfset(char*);
  
#define fsetperturbativeorder FC_FUNC(setperturbativeorder,SETPERTURBATIVEORDER)
  void fsetperturbativeorder(int*);

#define fgetperturbativeorder FC_FUNC(getperturbativeorder,GETPERTURBATIVEORDER)
  int fgetperturbativeorder();

#define fsetpolemasses FC_FUNC(setpolemasses,SETPOLEMASSES)
  void fsetpolemasses(double*,double*,double*);

#define fsetrenfacratio FC_FUNC(setrenfacratio,SETRENFACRATIO)
  void fsetrenfacratio(double*);

#define fsetreplica FC_FUNC(setreplica,SETREPLICA)
  void fsetreplica(int*);

#define fsettheory FC_FUNC(settheory, SETTHEORY)
  void fsettheory(char*);

#define fsetvfns FC_FUNC(setvfns,SETVFNS)
  void fsetvfns();

#define fgetapfelversion FC_FUNC(getapfelversion,GETAPFELVERSION)
  void fgetapfelversion(char*,int len);

#define fdisxsec FC_FUNC(fdisxsec,FDISXSEC)
  void fdisxsec(double*,double*,double*,double*,double*,char*,char*,int*,
                char*,int*,char*,char*,double*,double*,double*,double*);

#define finitializeapfel_dis FC_FUNC(initializeapfel_dis,INITIALIZEAPFEL_DIS)
  void finitializeapfel_dis(void);

#define fcomputestructurefunctionsapfel FC_FUNC(computestructurefunctionsapfel,COMPUTESTRUCTUREFUNCTIONSAPFEL)
  void fcomputestructurefunctionsapfel(double*, double*);

#define fsetmassscheme FC_FUNC(setmassscheme,SETMASSSCHEME)
  void fsetmassscheme(char*);

#define fsetpolarizationdis FC_FUNC(setpolarizationdis,SETPOLARIZATIONDIS)
  void fsetpolarizationdis(double*);

#define fsetprocessdis FC_FUNC(setprocessdis,SETPROCESSDIS)
  void fsetprocessdis(char*);

#define fsetprojectiledis FC_FUNC(setprojectiledis,SETPROJECTILEDIS)
  void fsetprojectiledis(char*);

#define fsettargetdis FC_FUNC(settargetdis,SETTARGETDIS)
  void fsettargetdis(char*);

#define fexternaldisoperator FC_FUNC(externaldisoperator,EXTERNALDISOPERATOR)
  double fexternaldisoperator(char*,int*,int*,double*,int*);

#define ff2light FC_FUNC(f2light,F2LIGHT)
  double ff2light(double*);

#define ff2charm FC_FUNC(f2charm,F2CHARM)
  double ff2charm(double*);

#define ff2bottom FC_FUNC(f2bottom,F2BOTTOM)
  double ff2bottom(double*);

#define ff2top FC_FUNC(f2top,F2TOP)
  double ff2top(double*);

#define ff2total FC_FUNC(f2total,F2TOTAL)
  double ff2total(double*);

#define ffllight FC_FUNC(fllight,FLLIGHT)
  double ffllight(double*);

#define fflcharm FC_FUNC(flcharm,FLCHARM)
  double fflcharm(double*);

#define fflbottom FC_FUNC(flbottom,FLBOTTOM)
  double fflbottom(double*);

#define ffltop FC_FUNC(fltop,FLTOP)
  double ffltop(double*);

#define ffltotal FC_FUNC(fltotal,FLTOTAL)
  double ffltotal(double*);

#define ff3light FC_FUNC(f3light,F3LIGHT)
  double ff3light(double*);

#define ff3charm FC_FUNC(f3charm,F3CHARM)
  double ff3charm(double*);

#define ff3bottom FC_FUNC(f3bottom,F3BOTTOM)
  double ff3bottom(double*);

#define ff3top FC_FUNC(f3top,F3TOP)
  double ff3top(double*);

#define ff3total FC_FUNC(f3total,F3TOTAL)
  double ff3total(double*);

#define fsetzmass FC_FUNC(setzmass,SETZMASS)
  void fsetzmass(double*);

#define fsetwmass FC_FUNC(setwmass,SETWMASS)
  void fsetwmass(double*);

#define fsetprotonmass FC_FUNC(setprotonmass,SETPROTONMASS)
  void fsetprotonmass(double*);

#define fsetsinthetaw FC_FUNC(setsinthetaw,SETSINTHETAW)
  void fsetsinthetaw(double*);

#define fsetckm FC_FUNC(setckm,SETCKM)
  void fsetckm(double*,double*,double*,double*,double*,double*,double*,double*,double*);

#define fsetgfermi FC_FUNC(setgfermi,SETGFERMI)
  void fsetgfermi(double*);

#define fgetzmass FC_FUNC(getzmass,GETZMASS)
  double fgetzmass();

#define fgetwmass FC_FUNC(getwmass,GETWMASS)
  double fgetwmass();

#define fgetprotonmass FC_FUNC(getprotonmass,GETPROTONMASS)
  double fgetprotonmass();

#define fgetsinthetaw FC_FUNC(getsinthetaw,GETSINTHETAW)
  double fgetsinthetaw();

#define fgetckm FC_FUNC(getckm,GETCKM)
  double fgetckm(int*,int*);

#define fgetgfermi FC_FUNC(fgetgfermi,FGETGFERMI)
  double fgetgfermi();

#define fenabletargetmasscorrections FC_FUNC(enabletargetmasscorrections,ENABLETARGETMASSCORRECTIONS)
  void fenabletargetmasscorrections(int*);

#define ffksimulator FC_FUNC(fksimulator,FKSIMULATOR)
  double ffksimulator(double*,double*,double*,int*,int*);

#define fsetfkobservable FC_FUNC(setfkobservable,SETFKOBSERVABLE)
  void fsetfkobservable(char*);

#define fgetfkobservable FC_FUNC(getfkobservable,GETFKOBSERVABLE)
  void fgetfkobservable();

#define ffkobservables FC_FUNC(fkobservables,FKOBSERVABLES)
  double ffkobservables(double*,double*,double*);

#define flistfunctions FC_FUNC(listfunctions,LISTFUNCTIONS)
  void flistfunctions(void);

}

#endif
