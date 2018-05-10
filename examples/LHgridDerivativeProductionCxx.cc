#include <iostream>
#include <cmath>
#include "APFEL/APFEL.h"

int main()
{
  APFEL::SetPerturbativeOrder(0);
  APFEL::SetPDFSet("NNPDF23_nlo_as_0118");

  APFEL::LHAPDFgridDerivative(0,"NNPDF23_nlo_as_0118_derived");

  return 0;
}
