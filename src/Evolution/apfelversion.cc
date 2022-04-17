// Define APFEL version string conversion
// Implements fake symbols for OS compilation

#include "APFEL/FortranWrappers.h"
#include <string>
#include <cstring>

#define STR_EXPAND(top) #top
#define STR(tok) STR_EXPAND(tok)

extern "C" {
  #define fgetapfelversion FC_FUNC(getapfelversion, GETAPFELVERSION)
  void fgetapfelversion(char* fversion, int) {
    std::string version = STR(APFEL_VERSION);
    strncpy(fversion, version.c_str(),  version.length()+1);
    for (size_t i = strlen(fversion); i < (unsigned) version.length()+1; ++i) {
      fversion[i] = ' ';
    }
  }
}

// define external functions for OS compilation
#ifndef DARWIN
   void externalsetapfel_(double, double, double*)  { return; }
   void externalsetapfel1_(double, double, double*) { return; }
   void externalsetapfelrep_(double, double, int, double*)  { return; }
   void externalsetapfelrep1_(double, double, int, double*) { return; }
   void externalsetapfellept_(double, double, int, double*, double*) { return; }
   void pretabulatedpdfsrep_(int, int, int, double*) { return; }
#endif
