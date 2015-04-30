#ifndef APFELfw_H
#define APFELfw_H

// Declarations for the Fortran/C interface

#include "APFEL/FortranWrappers.h"

extern "C" {
  

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

#define fselectcharge FC_FUNC(selectcharge,SELECTCHARGE)
  void fselectcharge(char*);

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

#define fcomputefktables FC_FUNC(computefktables,COMPUTEFKTABLES)
  void fcomputefktables(char*,char*,double*,int*);    

#define fcomputehardcrosssectionsdy FC_FUNC(computehardcrosssectionsdy,COMPUTEHARDCROSSSECTIONDY)
  void fcomputehardcrosssectionsdy(char*,char*);    

}

#endif