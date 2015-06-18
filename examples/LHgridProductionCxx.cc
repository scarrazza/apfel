#include <iostream>
#include <cmath>
#include "APFEL/APFEL.h"

int main()
{
  double Qin = sqrt(2.0);
  APFEL::LHAPDFgrid(0,Qin,"ApfelPDFs");

  return 0;
}
