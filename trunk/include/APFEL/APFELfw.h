#ifndef APFELfw_H
#define APFELfw_H

// Declarations for the Fortran/C interface

#include "APFEL/FortranWrappers.h"

extern "C" {
  
#define finitializeapfel FC_FUNC(initializeapfel,INITIALIZEAPFEL)
  void finitializeapfel(void);

#define fevolveapfel FC_FUNC(evolveapfel,EVOLVEAPFEL)
  void fevolveapfel(double*,double*);

#define fxpdf FC_FUNC(xpdf,XPDF)
  double fxpdf(int*, double*);

#define fxpdfj FC_FUNC(xpdfj,XPDFJ)
  double fxpdfj(int*, double*);

#define fxgamma FC_FUNC(xgamma,XGAMMA)
  double fxgamma(double*);

#define fxgammaj FC_FUNC(xgammaj,XGAMMAJ)
  double fxgammaj(double*);

#define flhapdfgrid FC_FUNC(lhapdfgrid,LHAPDFGRID)
  void flhapdfgrid(int*, double*, char*);

#define fexternalevolutionoperator FC_FUNC(fexternalevolutionoperator,FEXTERNALEVOLUTIONOPERATOR)
  void fexternalevolutionoperator(double*,double*,int*,double*,double*);

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

#define fcleanup FC_FUNC(cleanup,CLEANUP)
  void fcleanup(void);

#define fenablewelcomemessage FC_FUNC(enablewelcomemessage,ENABLEWELCOMEMESSAGE)
  void fenablewelcomemessage(bool*);

#define fenableevolutionoperator FC_FUNC(enableevolutionoperator,ENABLEEVOLUTIONOPERATOR)
  void fenableevolutionoperator(bool*);

#define fheavyquarkmass FC_FUNC(heavyquarkmass,HEAVYQUARKMASS)
  double fheavyquarkmass(int*,double*);

#define fsetalphaqcdref FC_FUNC(setalphaqcdref,SETALPHAQCDREF)
  void fsetalphaqcdref(double*,double*);

#define fsetalphaqedref FC_FUNC(setalphaqedref,SETALPHAQEDREF)
  void fsetalphaqedref(double*,double*);

#define fsetqlimits FC_FUNC(setqlimits,SETQLIMITS)
  void fsetqlimits(double*,double*);

#define fsetffns FC_FUNC(setffns,SETFFNS)
  void fsetffns(int*);

#define fsetgridparameters FC_FUNC(setgridparameters,SETGRIDPARAMETERS)
  void fsetgridparameters(int*,int*,int*,double*);

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

#define fsetpolemasses FC_FUNC(setpolemasses, SETPOLEMASSES)
  void fsetpolemasses(double*,double*,double*);

#define fsetrenfacratio FC_FUNC(setrenfacratio, SETRENFACRATIO)
  void fsetrenfacratio(double*);

#define fsetreplica FC_FUNC(setreplica, SETREPLICA)
  void fsetreplica(int*);

#define fsettheory FC_FUNC(settheory, SETTHEORY)
  void fsettheory(char*);

#define fsetvfns FC_FUNC(setvfns, SETVFNS)
  void fsetvfns();

#define fgetapfelversion FC_FUNC(getapfelversion, GETAPFELVERSION)
  void fgetapfelversion(char*,int len);

#define fdisxsec FC_FUNC(fdisxsec, FDISXSEC)
  void fdisxsec(double*,double*,double*,double*,double*,char*,char*,int*,
                char*,int*,char*,char*,double*,double*,double*,double*);
}

#endif
