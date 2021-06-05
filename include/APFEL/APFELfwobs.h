#ifndef APFELfw_H
#define APFELfw_H

// Declarations for the Fortran/C interface

#include "APFEL/FortranWrappers.h"

extern "C" {
  
#define finitializeapfel_dis FC_FUNC(initializeapfel_dis,INITIALIZEAPFEL_DIS)
  void finitializeapfel_dis(void);

#define fcomputestructurefunctionsapfel FC_FUNC(computestructurefunctionsapfel,COMPUTESTRUCTUREFUNCTIONSAPFEL)
  void fcomputestructurefunctionsapfel(double*, double*);

#define fcachestructurefunctionsapfel FC_FUNC(cachestructurefunctionsapfel,CACHESTRUCTUREFUNCTIONSAPFEL)
  void fcachestructurefunctionsapfel(double*);

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

#define fg1light FC_FUNC(g1light,G1LIGHT)
  double fg1light(double*);

#define fg1charm FC_FUNC(g1charm,G1CHARM)
  double fg1charm(double*);

#define fg1bottom FC_FUNC(g1bottom,G1BOTTOM)
  double fg1bottom(double*);

#define fg1top FC_FUNC(g1top,G1TOP)
  double fg1top(double*);

#define fg1total FC_FUNC(g1total,G1TOTAL)
  double fg1total(double*);

#define fgllight FC_FUNC(gllight,GLLIGHT)
  double fgllight(double*);

#define fglcharm FC_FUNC(glcharm,GLCHARM)
  double fglcharm(double*);

#define fglbottom FC_FUNC(glbottom,GLBOTTOM)
  double fglbottom(double*);

#define fgltop FC_FUNC(gltop,GLTOP)
  double fgltop(double*);

#define fgltotal FC_FUNC(gltotal,GLTOTAL)
  double fgltotal(double*);

#define fg4light FC_FUNC(g4light,G4LIGHT)
  double fg4light(double*);

#define fg4charm FC_FUNC(g4charm,G4CHARM)
  double fg4charm(double*);

#define fg4bottom FC_FUNC(g4bottom,G4BOTTOM)
  double fg4bottom(double*);

#define fg4top FC_FUNC(g4top,G4TOP)
  double fg4top(double*);

#define fg4total FC_FUNC(g4total,G4TOTAL)
  double fg4total(double*);

#define fstructurefunctionxq FC_FUNC(structurefunctionxq,STRUCTUREFUNCTIONXQ)
  double fstructurefunctionxq(char*,char*,char*,double*,double*);

#define fsetzmass FC_FUNC(setzmass,SETZMASS)
  void fsetzmass(double*);

#define fsetwmass FC_FUNC(setwmass,SETWMASS)
  void fsetwmass(double*);

#define fsetprotonmass FC_FUNC(setprotonmass,SETPROTONMASS)
  void fsetprotonmass(double*);

#define fsetsin2thetaw FC_FUNC(setsin2thetaw,SETSIN2THETAW)
  void fsetsin2thetaw(double*);

#define fsetckm FC_FUNC(setckm,SETCKM)
  void fsetckm(double*,double*,double*,double*,double*,double*,double*,double*,double*);

#define fsetpropagatorcorrection FC_FUNC(setpropagatorcorrection,SETPROPAGATORCORRECTION)
  void fsetpropagatorcorrection(double*);

#define fsetewcouplings FC_FUNC(setewcouplings,SETEWCOUPLINGS)
  void fsetewcouplings(double*,double*,double*,double*);

#define fsetgfermi FC_FUNC(setgfermi,SETGFERMI)
  void fsetgfermi(double*);

#define fsetrenqratio FC_FUNC(setrenqratio,SETRENQRATIO)
  void fsetrenqratio(double*);

#define fsetfacqratio FC_FUNC(setfacqratio,SETFACQRATIO)
  void fsetfacqratio(double*);

#define fenabledynamicalscalevariations FC_FUNC(enabledynamicalscalevariations,ENABLEDYNAMICALSCALEVARIATIONS)
  void fenabledynamicalscalevariations(int*);

#define fenableintrinsiccharm FC_FUNC(enableintrinsiccharm,ENABLEINTRINSICCHARM)
  void fenableintrinsiccharm(int*);

#define fgetzmass FC_FUNC(getzmass,GETZMASS)
  double fgetzmass();

#define fgetwmass FC_FUNC(getwmass,GETWMASS)
  double fgetwmass();

#define fgetprotonmass FC_FUNC(getprotonmass,GETPROTONMASS)
  double fgetprotonmass();

#define fgetsin2thetaw FC_FUNC(getsin2thetaw,GETSIN2THETAW)
  double fgetsin2thetaw();

#define fgetckm FC_FUNC(getckm,GETCKM)
  double fgetckm(int*,int*);

#define fgetgfermi FC_FUNC(getgfermi,GETGFERMI)
  double fgetgfermi();

#define fgetsiatotalcrosssection FC_FUNC(getsiatotalcrosssection,GETSIATOTALCROSSSECTION)
  double fgetsiatotalcrosssection(int*,double*,char*);

#define fenabletargetmasscorrections FC_FUNC(enabletargetmasscorrections,ENABLETARGETMASSCORRECTIONS)
  void fenabletargetmasscorrections(int*);

#define fenabledampingfonll FC_FUNC(enabledampingfonll,ENABLEDAMPINGFONLL)
  void fenabledampingfonll(int*);

#define fsetdampingpowerfonll FC_FUNC(setdampingpowerfonll,SETDAMPINGPOWERFONLL)
  void fsetdampingpowerfonll(int*,int*,int*);

#define fcomputechargesdis FC_FUNC(computechargesdis,COMPUTECHARGESDIS)
  void fcomputechargesdis(double*,double*,double*,double*);

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

#define fenablesfnloqedcorrections FC_FUNC(enablesfnloqedcorrections,ENABLESFNLOQEDCORRECTIONS)
  void fenablesfnloqedcorrections(int*);

#define flhapdfgridstructurefunctions FC_FUNC(lhapdfgridstructurefunctions,LHAPDFGRIDSTRUCTUREFUNCTIONS)
  void flhapdfgridstructurefunctions(int*, double*, char*);

#define fsetscalevariationprocedure FC_FUNC(setscalevariationprocedure,SETSCALEVARIATIONPROCEDURE)
  void fsetscalevariationprocedure(int*);
}

#endif
