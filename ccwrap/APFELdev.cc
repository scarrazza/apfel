// C++ definition

#include "APFEL/APFELdev.h"
#include "APFEL/APFELevol.h"
#include "APFEL/APFELobs.h"
#include <cstdlib>
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
    if (!atoi(par.at(kQED).c_str()))
      APFEL::SetTheory(string("QCD"));
    else
      {
	APFEL::SetTheory(string("QUniD"));
	APFEL::EnableNLOQEDCorrections(true);
      }
    APFEL::SetPerturbativeOrder(atoi(par.at(kPTO).c_str()));

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
    APFEL::SetAlphaQCDRef(atof(par.at(kalphas).c_str()), atof(par.at(kQref).c_str()));
    if (atoi(par.at(kQED).c_str())) APFEL::SetAlphaQEDRef(atof(par.at(kalphaqed).c_str()),atof(par.at(kQedref).c_str()));
    
    // EW
    APFEL::SetWMass(atof(par.at(kMW).c_str()));
    APFEL::SetZMass(atof(par.at(kMZ).c_str()));
    APFEL::SetGFermi(atof(par.at(kGF).c_str()));
    APFEL::SetSin2ThetaW(atof(par.at(kSIN2TW).c_str()));

    vector<string> ckm = split(par.at(kCKM));

    APFEL::SetCKM(atof(ckm[0].c_str()), atof(ckm[1].c_str()), atof(ckm[2].c_str()),
                  atof(ckm[3].c_str()), atof(ckm[4].c_str()), atof(ckm[5].c_str()),
                  atof(ckm[6].c_str()), atof(ckm[7].c_str()), atof(ckm[8].c_str()));

    // TMCs
    APFEL::SetProtonMass(atof(par.at(kMP).c_str()));
    if (atoi(par.at(kTMC).c_str())) APFEL::EnableTargetMassCorrections(true);

    // Heavy Quark Masses
    if (par.at(kHQ).compare("POLE") == 0 )
      APFEL::SetPoleMasses(atof(par.at(kmc).c_str()), atof(par.at(kmb).c_str()), atof(par.at(kmt).c_str()));
    else if (par.at(kHQ).compare("MSBAR") == 0 )
    {
      APFEL::SetMSbarMasses(atof(par.at(kmc).c_str()), atof(par.at(kmb).c_str()), atof(par.at(kmt).c_str()));
      APFEL::SetMassScaleReference(atof(par.at(kQmc).c_str()), atof(par.at(kQmb).c_str()), atof(par.at(kQmt).c_str()));
    }
    else
    {
      cerr << "Error: Unrecognised HQMASS"<<endl;
      exit(-1);
    }

    // Heavy Quark schemes
    APFEL::SetMassScheme(par.at(kFNS));
    APFEL::EnableDampingFONLL(atoi(par.at(kDAMP).c_str()));
    if (par.at(kFNS).compare("FFNS") == 0)
      APFEL::SetFFNS(atoi(par.at(kNfFF).c_str()));
    else
      APFEL::SetVFNS();
    
    APFEL::SetMaxFlavourAlpha(atoi(par.at(kMaxNfAs).c_str()));
    APFEL::SetMaxFlavourPDFs(atoi(par.at(kMaxNfPdf).c_str()));

    // Scale ratios
    APFEL::SetRenFacRatio(atof(par.at(kXIR).c_str())/atof(par.at(kXIF).c_str()));
    APFEL::SetRenQRatio(atof(par.at(kXIR).c_str()));
    APFEL::SetFacQRatio(atof(par.at(kXIF).c_str()));

    // Small-x resummation
    APFEL::SetSmallxResummation(atoi(par.at(kSxRes).c_str()), par.at(kSxOrd));
    APFEL::SetMassMatchingScales(atof(par.at(kcThr).c_str()),atof(par.at(kbThr).c_str()),atof(par.at(ktThr).c_str()));

    // Intrinsic charm
    APFEL::EnableIntrinsicCharm(atoi(par.at(kIC).c_str()));

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
