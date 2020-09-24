%module apfel

%include "std_string.i"
%include "std_vector.i"
%include "carrays.i"
%array_functions(double,doubles)
//%include "std_except.i"

%{
  #define SWIG_FILE_WITH_INIT
  #include "APFEL/APFELevol.h"
  #include "APFEL/APFELobs.h"
  #include <cstddef>
%}

namespace APFEL {}

%include "APFEL/APFELevol.h"
%include "APFEL/APFELobs.h"
