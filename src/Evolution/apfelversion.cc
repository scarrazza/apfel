#include "APFEL/FortranWrappers.h"
#include <string>
#include <cstring>
#include <cassert>

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <sys/stat.h>

#ifndef NOLHAPDF
#include "LHAPDF/LHAPDF.h"
#endif

using namespace std;

#define STR_EXPAND(top) #top
#define STR(tok) STR_EXPAND(tok)

extern "C" {

  #define fgetapfelversion FC_FUNC(getapfelversion, GETAPFELVERSION)
  void fgetapfelversion(char* fversion, int length) {
    string version = STR(APFEL_VERSION);    
    strncpy(fversion, version.c_str(),  version.length()+1);    
    for (size_t i = strlen(fversion); i < (unsigned) version.length()+1; ++i) {
      fversion[i] = ' ';
    }
  }

  bool islhapdf6_()
  {
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
    return true;
#else
    return false;
#endif
  }

  void mkdir_(char* name, int length)
  {
    mkdir(name, 0777);
  }

#ifdef NOLHAPDF
  void stop()
  {
    printf(" [Error] LHAPDF support disabled, please recompile APFEL without --disable-lhapdf.\n");
    exit(-1);
  }

  void getqmass_(){ stop(); }
  void initpdf_(){ stop(); }
  void alphaspdf_(){ stop(); }
  void has_photon_(){ stop(); }
  void initpdfsetbyname_(){ stop(); }
  void evolvepdfphoton_(){ stop(); }
  void evolvepdf_(){ stop(); }
  void numberpdf_(){ stop(); }
#endif

  // remove if mac
#ifndef DARWIN
  void externalsetapfel_(double x, double Q, double* xf)
  {
    return;
  }

  void externalsetapfel1_(double x, double Q, double* xf)
  {
    return;
  }

  void externalsetapfellept_(double x, double Q, int irep, 
			     double* xl, double* xf)
  {
    return;
  }
#endif

}
