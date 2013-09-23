%module apfel

%include "std_string.i"
%include "std_vector.i"
//%include "std_except.i"

%{
  #define SWIG_FILE_WITH_INIT
  #include "APFEL/APFEL.h"
  #include <cstddef>
%}

namespace APFEL {}

%include "APFEL/APFEL.h"
