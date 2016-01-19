#ifndef APFELINIT_H
#define APFELINIT_H

#include <map>
#include <vector>
#include <string>
#include <iterator>
using std::map;
using std::vector;
using std::string;

namespace APFEL
{
  /**
   * @brief A container for the APFEL setup parameters
   *
   * Data structure for one time initialization for developers only.
   */
  const string kID       = "ID";
  const string kPTO      = "PTO";
  const string kFNS      = "FNS";
  const string kDAMP     = "DAMP";
  const string kIC       = "IC";
  const string kModEv    = "ModEv";
  const string kXIR      = "XIR";
  const string kXIF      = "XIF";
  const string kNfFF     = "NfFF";
  const string kMaxNfAs  = "MaxNfAs";
  const string kMaxNfPdf = "MaxNfPdf";
  const string kQ0       = "Q0";
  const string kalphas   = "alphas";
  const string kQref     = "Qref";
  const string kQED      = "QED";
  const string kalphaqed = "alphaqed";
  const string kQedref   = "Qedref";
  const string kSxRes    = "SxRes";
  const string kSxOrd    = "SxOrd";
  const string kHQ       = "HQ";
  const string kmc       = "mc";
  const string kQmc      = "Qmc";
  const string kmb       = "mb";
  const string kQmb      = "Qmb";
  const string kmt       = "mt";
  const string kQmt      = "Qmt";
  const string kcThr     = "kcThr";
  const string kbThr     = "kbThr";
  const string ktThr     = "ktThr";
  const string kCKM      = "CKM";
  const string kMZ       = "MZ";
  const string kMW       = "MW";
  const string kGF       = "GF";
  const string kSIN2TW   = "SIN2TW";
  const string kTMC      = "TMC";
  const string kMP       = "MP";
  const string kComments = "Comments";

  const string values[] = { kID, kPTO, kFNS, kDAMP, kIC, kModEv, kXIR, kXIF, kNfFF, kMaxNfAs,
                            kMaxNfPdf, kQ0, kalphas, kQref, kQED, kalphaqed, kQedref, kSxRes,
                            kSxOrd, kHQ, kmc, kQmc, kmb, kQmb, kmt, kQmt, kcThr, kbThr, ktThr,
			    kCKM, kMZ, kMW, kGF, kSIN2TW, kTMC, kMP, kComments };

  const vector<string> kValues(values,  values + sizeof values/sizeof values[0]);

  //! Setup APFEL with the previous data structure - first item is the string key see kvalues
  void SetParam(map<string,string> const& par);
}

#endif
