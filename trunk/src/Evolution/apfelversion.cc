#include "APFEL/FortranWrappers.h"
#include <string>
#include <cstring>
#include <cassert>

#include <iostream>

using namespace std;

#ifndef APFEL_VERSION
#define APFEL_VERSION "2.0.2"
#endif


extern "C" {

  #define fgetapfelversion FC_FUNC(getapfelversion, GETAPFELVERSION)
  void fgetapfelversion(char* fversion, int length) {
    string version = APFEL_VERSION;    
    strncpy(fversion, version.c_str(),  version.length()+1);    
    for (size_t i = strlen(fversion); i < (unsigned) version.length()+1; ++i) {
      fversion[i] = ' ';
    }
  }

}
