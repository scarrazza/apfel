#include "APFEL/FortranWrappers.h"
#include <string>
#include <cstring>
#include <cassert>

using namespace std;

#ifndef APFEL_VERSION
#define APFEL_VERSION "1.0.1"
#endif


extern "C" {

  #define fgetapfelversion FC_FUNC(getapfelversion, GETAPFELVERSION)
  void fgetapfelversion(char* fversion, int length) {
    string version = APFEL_VERSION;
    strncpy(fversion, version.c_str(), length);
    for (size_t i = strlen(fversion); i < (unsigned) length; ++i) {
      fversion[i] = ' ';
    }
  }

}
