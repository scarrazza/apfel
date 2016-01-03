// C++ definition

#include "APFEL/APFELinit.h"
#include "APFEL/APFELevol.h"
#include "APFEL/APFELobs.h"
#include <iostream>
using std::cerr;
using std::endl;

namespace APFEL
{
  void SetParam(apfel_param const& par)
  {
    // Cleanup APFEL common blocks
    APFEL::CleanUp();

    // Theory, perturbative order of evolution
    if (!par.QED)
      APFEL::SetTheory(string("QCD"));
    else
      APFEL::SetTheory(string("QUniD"));
    APFEL::SetPerturbativeOrder(par.evol_pto);

    if (par.MODEV.compare("EXA") == 0)
    {
      APFEL::SetPDFEvolution("exactalpha");
      APFEL::SetAlphaEvolution("exact");
    }
    else if (par.MODEV.compare("EXP") == 0)
    {
      APFEL::SetPDFEvolution("expandalpha");
      APFEL::SetAlphaEvolution("expanded");
    }
    else if (par.MODEV.compare("TRN") == 0)
    {
      APFEL::SetPDFEvolution("truncated");
      APFEL::SetAlphaEvolution("expanded");
    }
    else
    {
      std::cerr << " ERROR: Unrecognised MODEV: "<<par.MODEV<<std::endl;
      exit(-1);
    }
    
    // Coupling
    APFEL::SetAlphaQCDRef(par.alphas, par.QREF);
    if (par.QED) APFEL::SetAlphaQEDRef(par.alpha_qed,par.QEDREF);
    
    // EW
    APFEL::SetWMass(par.mw);
    APFEL::SetZMass(par.mz);
    APFEL::SetGFermi(par.gf);
    APFEL::SetSin2ThetaW(par.sin2tw);

    APFEL::SetCKM(par.CKM[0][0], par.CKM[0][1], par.CKM[0][2],
                  par.CKM[1][0], par.CKM[1][1], par.CKM[1][2],
                  par.CKM[2][0], par.CKM[2][1], par.CKM[2][2]);

    // TMCs
    APFEL::SetProtonMass(par.Mt);
    if (par.TMC) APFEL::EnableTargetMassCorrections(true);

    // Heavy Quark Masses
    if (par.HQMASS.compare("POLE") == 0 ) 
      APFEL::SetPoleMasses(par.mc, par.mb, par.mt);
    else if (par.HQMASS.compare("MSBAR") == 0 )
    {
      APFEL::SetMSbarMasses(par.mc, par.mb, par.mt);
      APFEL::SetMassScaleReference(par.Qmc, par.Qmb, par.Qmt);
    }
    else
    {
      cerr << "Error: Unrecognised HQMASS"<<endl;
      exit(-1);
    }
    
    // Heavy Quark schemes
    APFEL::SetMassScheme(par.FNS);
    APFEL::EnableDampingFONLL(par.damp);
    if (par.FNS.compare("FFNS") == 0)
      APFEL::SetFFNS(par.nff);
    else
      APFEL::SetVFNS();
    
    APFEL::SetMaxFlavourAlpha(par.nf_as);
    APFEL::SetMaxFlavourPDFs(par.nf_pdf);
    
    // Truncated Epsilon
    APFEL::SetEpsilonTruncation(par.truc_eps);

    // Scale ratios
    APFEL::SetRenFacRatio(par.xiR/par.xiF);
    APFEL::SetRenQRatio(par.xiR);
    APFEL::SetFacQRatio(par.xiF);

    // Small-x resummation
    APFEL::SetSmallxResummation(par.SxRes, par.SxOrd);
    
    // Set maximum scale
    APFEL::SetQLimits(par.Q0, par.QM );

    if (par.SIA) 
    {
      APFEL::SetPDFSet("kretzer");
      APFEL::SetTimeLikeEvolution(true);
    }

    return;
  };
}; 
