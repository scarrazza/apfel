#ifndef COMMON_H
#define COMMON_H

#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TString.h"
#include "TLegend.h"

#include <QString>

const double lharanges[][6] = {
  {1e-5, 1.0, -0.1, 1.3, -0.5, 0.5}, // tbar
  {1e-5, 1.0, -0.1, 1.3, -0.5, 0.5}, // bbar
  {1e-5, 1.0, -1.0, 1.0, -0.1, 0.1}, // cbar
  {1e-5, 1.0, -1.5, 2.5, -0.1, 0.1}, // sbar
  {1e-5, 1.0, -0.1, 1.3, -0.1, 0.4}, // ubar
  {1e-5, 1.0, -0.1, 1.3, -0.1, 0.4}, // dbar
  {1e-5, 1.0, -2.5, 7.0, -0.1, 0.6}, // g
  {1e-5, 1.0, -0.1, 1.3, -0.1, 0.6}, // d
  {1e-5, 1.0, -0.1, 1.3, -0.1, 1.0}, // u
  {1e-5, 1.0, -1.5, 2.5, -0.1, 0.1}, // s
  {1e-3, 1.0, -0.1, 1.0, -0.1, 0.1}, // c
  {1e-3, 1.0, -0.1, 1.0, -0.1, 0.5}, // b
  {1e-5, 1.0, -0.1, 1.3, -0.5, 0.5}, // t
  {1e-5, 1.0, -3.0, 3.0, -1.0, 1.0}, // photon
                        };

const QString name[] = {
  "x#bar{t}(x,Q)",
  "x#bar{b}(x,Q)",
  "x#bar{c}(x,Q)",
  "x#bar{s}(x,Q)",
  "x#bar{u}(x,Q)",
  "x#bar{d}(x,Q)",
  "xg(x,Q)",
  "xd(x,Q)",
  "xu(x,Q)",
  "xs(x,Q)",
  "xc(x,Q)",
  "xb(x,Q)",
  "xt(x,Q)",
  "x#gamma(x,Q)"
};

const int colors[] = { kGreen, kBlue, kRed, kOrange,
                 kViolet, kMagenta, kBlack, kYellow,
                 kCyan, kGray};

double ComputeAVG(int n, int ix, double **x);

double ComputeStdDev(int n, int ix, double **x);


#endif // COMMON_H
