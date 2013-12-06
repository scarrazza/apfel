#ifndef COMMON_H
#define COMMON_H

#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
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

const int colors2[] = {kYellow, kBlack, kMagenta, kCyan+1, kBlue+1, kGreen+1,
                       kRed,
                       kGreen, kBlue, kCyan+1, kMagenta, kBlack, kYellow};

const int fillcolor[] = {kGreen, kBlue, kRed, kCyan, kViolet,kYellow,kOrange,kAzure,kBlack,kGray,
                        kGreen+1, kBlue+1, kRed+1, kCyan+1, kViolet+1,kYellow+1,kOrange+1,kAzure+1,kBlack+1,kGray+1,
                        kGreen+2, kBlue+2, kRed+2, kCyan+2, kViolet+2,kYellow+2,kOrange+2,kAzure+2,kBlack+2,kGray+2,
                        kGreen+3, kBlue+3, kRed+3, kCyan+3, kViolet+3,kYellow+3,kOrange+3,kAzure+3,kBlack+3,kGray+3};

const int fillStyle[] = {1001, 3005, 3004, 3006, 3007, 3001, 3002, 3003, 3004, 3005,
                         3005, 3004, 3006, 3007, 3001, 3002, 3003, 3004, 3005,
                         3005, 3004, 3006, 3007, 3001, 3002, 3003, 3004, 3005,
                         3005, 3004, 3006, 3007, 3001, 3002, 3003, 3004, 3005};

#endif // COMMON_H
