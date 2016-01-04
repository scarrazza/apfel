// C++ definition

#include "APFEL/APFELdev.h"
#include "APFEL/APFELevol.h"
#include "APFEL/APFELobs.h"
#include <iostream>
#include <sstream>
using std::cerr;
using std::endl;
using std::cout;
using std::stringstream;
using std::istream_iterator;

namespace APFEL
{
  vector<string> split(string const &input)
  {
        stringstream strstr(input);
        istream_iterator<string> it(strstr);
        istream_iterator<string> end;
        vector<string> results(it, end);
        return results;
  }

  void SetParam(map<string,string> const& par)
  {
    // Cleanup APFEL common blocks
    APFEL::CleanUp();

    // Theory, perturbative order of evolution
    if (!stoi(par.at(kQED)))
      APFEL::SetTheory(string("QCD"));
    else
      APFEL::SetTheory(string("QUniD"));
    APFEL::SetPerturbativeOrder(stoi(par.at(kPTO)));

    if (par.at(kModEv).compare("EXA") == 0)
    {
      APFEL::SetPDFEvolution("exactalpha");
      APFEL::SetAlphaEvolution("exact");
    }
    else if (par.at(kModEv).compare("EXP") == 0)
    {
      APFEL::SetPDFEvolution("expandalpha");
      APFEL::SetAlphaEvolution("expanded");
    }
    else if (par.at(kModEv).compare("TRN") == 0)
    {
      APFEL::SetPDFEvolution("truncated");
      APFEL::SetAlphaEvolution("expanded");
    }
    else
    {
      std::cerr << " ERROR: Unrecognised MODEV: "<< par.at(kModEv) <<std::endl;
      exit(-1);
    }
        
    // Coupling
    APFEL::SetAlphaQCDRef(stod(par.at(kalphas)), stod(par.at(kQref)));
    if (stoi(par.at(kQED))) APFEL::SetAlphaQEDRef(stod(par.at(kalphaqed)),stod(par.at(kQedref)));
    
    // EW
    APFEL::SetWMass(stod(par.at(kMW)));
    APFEL::SetZMass(stod(par.at(kMZ)));
    APFEL::SetGFermi(stod(par.at(kGF)));
    APFEL::SetSin2ThetaW(stod(par.at(kSIN2TW)));

    vector<string> ckm = split(par.at(kCKM));

    APFEL::SetCKM(stod(ckm[0]), stod(ckm[1]), stod(ckm[2]),
                  stod(ckm[3]), stod(ckm[4]), stod(ckm[5]),
                  stod(ckm[6]), stod(ckm[7]), stod(ckm[8]));

    // TMCs
    APFEL::SetProtonMass(stod(par.at(kMP)));
    if (stoi(par.at(kTMC))) APFEL::EnableTargetMassCorrections(true);

    // Heavy Quark Masses
    if (par.at(kHQ).compare("POLE") == 0 )
      APFEL::SetPoleMasses(stod(par.at(kmc)), stod(par.at(kmb)), stod(par.at(kmt)));
    else if (par.at(kHQ).compare("MSBAR") == 0 )
    {
      APFEL::SetMSbarMasses(stod(par.at(kmc)), stod(par.at(kmb)), stod(par.at(kmt)));
      APFEL::SetMassScaleReference(stod(par.at(kQmc)), stod(par.at(kQmb)), stod(par.at(kQmt)));
    }
    else
    {
      cerr << "Error: Unrecognised HQMASS"<<endl;
      exit(-1);
    }

    // Heavy Quark schemes
    APFEL::SetMassScheme(par.at(kFNS));
    APFEL::EnableDampingFONLL(stoi(par.at(kDAMP)));
    if (par.at(kFNS).compare("FFNS") == 0)
      APFEL::SetFFNS(stoi(par.at(kNfFF)));
    else
      APFEL::SetVFNS();
    
    APFEL::SetMaxFlavourAlpha(stoi(par.at(kMaxNfAs)));
    APFEL::SetMaxFlavourPDFs(stoi(par.at(kMaxNfPdf)));

    // Scale ratios
    APFEL::SetRenFacRatio(stod(par.at(kXIR))/stod(par.at(kXIF)));
    APFEL::SetRenQRatio(stod(par.at(kXIR)));
    APFEL::SetFacQRatio(stod(par.at(kXIF)));

    // Small-x resummation
    APFEL::SetSmallxResummation(stoi(par.at(kSxRes)), par.at(kSxOrd));

    // Not included in the map
    /*
    // Truncated Epsilon
    APFEL::SetEpsilonTruncation(1E-1);

    // Set maximum scale
    APFEL::SetQLimits(par.Q0, par.QM );

    if (par.SIA)
    {
      APFEL::SetPDFSet("kretzer");
      APFEL::SetTimeLikeEvolution(true);
    }
    */

    return;
  }
}
