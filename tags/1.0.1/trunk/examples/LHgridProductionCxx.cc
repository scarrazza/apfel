#include <iostream>
#include <cmath>
#include "APFEL/APFEL.h"

int main()
{
  APFEL::SetPoleMasses(sqrt(2.0),1e5,1e5);
  double Qin = sqrt(2.0);

  APFEL::LHAPDFgrid(0,Qin,"ApfelPDFs");

  return 0;
}
