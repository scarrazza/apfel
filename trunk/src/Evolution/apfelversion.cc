#include "APFEL/FortranWrappers.h"
#include "LHAPDF/LHAPDF.h"
#include <string>
#include <cstring>
#include <cassert>

#include <iostream>
#include <cstdlib>
#include <sys/stat.h>

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

}
