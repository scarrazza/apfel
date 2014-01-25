#ifndef APFELfw_H
#define APFELfw_H

// Declarations for the Fortran/C interface

#include "APFEL/FortranWrappers.h"

extern "C" {
  
#define finitializeapfel FC_FUNC(finitializeapfel, FINITIALIZEAPFEL)
  void finitializeapfel(void);

#define fevolveapfel FC_FUNC(fevolveapfel, FEVOLVEAPFEL)
  void fevolveapfel(double*,double*);

#define fxpdf FC_FUNC(fxpdf, FXPDF)
  double fxpdf(int*, double*);

#define fxpdfj FC_FUNC(fxpdfj, FXPDFJ)
  double fxpdfj(int*, double*);

#define fxgamma FC_FUNC(fxgamma, FXGAMMA)
  double fxgamma(double*);

#define fxgammaj FC_FUNC(fxgammaj, FXGAMMAJ)
  double fxgammaj(double*);

#define flhapdfgrid FC_FUNC(flhapdfgrid, FLHAPDFGRID)
  void flhapdfgrid(int*, double*, char*);

#define fexternalevolutionoperator FC_FUNC(fexternalevolutionoperator, FEXTERNALEVOLUTIONOPERATOR)
  void fexternalevolutionoperator(double*,double*,int*,double*,double*);

#define falphaqcd FC_FUNC(falphaqcd, FALPHAQCD)
  double falphaqcd(double*);

#define falphaqed FC_FUNC(falphaqed, FALPHAQED)
  double falphaqed(double*);

#define fnpdf FC_FUNC(fnpdf,FNPDF)
  double fnpdf(int*,int*);

#define fngamma FC_FUNC(fngamma,FNGAMMA)
  double fngamma(int*);

#define flumi FC_FUNC(flumi,FLUMI)
  double flumi(int*,int*,double*);

#define fcleanup FC_FUNC(fcleanup,FCLEANUP)
  void fcleanup(void);

#define fenablewelcomemessage FC_FUNC(fenablewelcomemessage,FENABLEWELCOMEMESSAGE)
  void fenablewelcomemessage(bool*);

#define fenableevolutionoperator FC_FUNC(fenableevolutionoperator,FENABLEEVOLUTIONOPERATOR)
  void fenableevolutionoperator(bool*);

#define fheavyquarkmass FC_FUNC(fheavyquarkmass,FHEAVYQUARKMASS)
  double fheavyquarkmass(int*,double*);

#define fsetalphaqcdref FC_FUNC(fsetalphaqcdref, FSETALPHAQCDREF)
  void fsetalphaqcdref(double*,double*);

#define fsetalphaqedref FC_FUNC(fsetalphaqedref, FSETALPHAQEDREF)
  void fsetalphaqedref(double*,double*);

#define fsetqlimits FC_FUNC(fsetqlimits, FSETQLIMITS)
  void fsetqlimits(double*,double*);

#define fsetffns FC_FUNC(fsetffns, FSETFFNS)
  void fsetffns(int*);

#define fsetgridparameters FC_FUNC(fsetgridparameters, FSETGRIDPARAMETERS)
  void fsetgridparameters(int*,int*,int*,double*);

#define fsetmaxflavouralpha FC_FUNC(fsetmaxflavouralpha, FSETMAXFLAVOURALPHA)
  void fsetmaxflavouralpha(int*);

#define fsetmaxflavourpdfs FC_FUNC(fsetmaxflavourpdfs, FSETMAXFLAVOURPDFS)
  void fsetmaxflavourpdfs(int*);

#define fsetmsbarmasses FC_FUNC(fsetmsbarmasses, FSETMSBARMASSES)
  void fsetmsbarmasses(double*,double*,double*);

#define fsetnumberofgrids FC_FUNC(fsetnumberofgrids, FSETNUMBEROFGRIDS)
  void fsetnumberofgrids(int*);

#define fsetpdfset FC_FUNC(fsetpdfset, FSETPDFSET)
  void fsetpdfset(char*);
  
#define fsetperturbativeorder FC_FUNC(fsetperturbativeorder,FSETPERTURBATIVEORDER)
  void fsetperturbativeorder(int*);

#define fsetpolemasses FC_FUNC(fsetpolemasses, FSETPOLEMASSES)
  void fsetpolemasses(double*,double*,double*);

#define fsetrenfacratio FC_FUNC(fsetrenfacratio, FSETRENFACRATIO)
  void fsetrenfacratio(double*);

#define fsetreplica FC_FUNC(fsetreplica, FSETREPLICA)
  void fsetreplica(int*);

#define fsettheory FC_FUNC(fsettheory, FSETTHEORY)
  void fsettheory(char*);

#define fsetvfns FC_FUNC(fsetvfns, FSETVFNS)
  void fsetvfns();

#define fgetapfelversion FC_FUNC(fgetapfelversion, FGETAPFELVERSION)
  void fgetapfelversion(char*,int len);

#define fdisxsec FC_FUNC(fdisxsec, FDISXSEC)
  void fdisxsec(double*,double*,double*,double*,double*,char*,char*,int*,
                char*,int*,char*,char*,double*,double*,double*,double*);
}

#endif
