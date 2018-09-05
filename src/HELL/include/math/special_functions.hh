#pragma once

namespace HELLx {
  const double ZETA2 = 1.6449340668482264;
  const double ZETA3 = 1.2020569031595942855;
  const double ZETA4 = 1.082323233711138191516;
  extern double Li2(double x);
  extern double Li3(double x);
  extern double Li4(double x);
  extern double HPLmp(double x);
  extern double ArcCoth(double x);
  extern double ArcCsch(double x);
  extern double dpsi(double Z, int M=1);
};

