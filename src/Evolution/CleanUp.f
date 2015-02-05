************************************************************************
*
*     CleanUp.f:
*
*     It unsets all the evolution parameters so that they need to be 
*     reset.
*
************************************************************************
      subroutine CleanUp
*
      implicit none
*
      include "../commons/Welcome.h"
      include "../commons/scales.h"
      include "../commons/Evs.h"
      include "../commons/Nf_FF.h"
      include "../commons/ipt.h"
      include "../commons/Th.h"
      include "../commons/alpha_ref_QCD.h"
      include "../commons/alpha_ref_QED.h"
      include "../commons/kren.h"
      include "../commons/mass_scheme.h"
      include "../commons/m2th.h"
      include "../commons/MaxFlavourPDFs.h"
      include "../commons/MaxFlavourAlpha.h"
      include "../commons/grid.h"
      include "../commons/pdfset.h"
      include "../commons/Replica.h"
      include "../commons/lock.h"
      include "../commons/EvolOp.h"
      include "../commons/lambda_ref_QCD.h"
      include "../commons/AlphaEvolution.h"
      include "../commons/PDFEvolution.h"
      include "../commons/TimeLike.h"
      include "../commons/Smallx.h"
      include "../commons/FastEvol.h"
      include "../commons/MassScheme.h"
      include "../commons/ProcessDIS.h"
      include "../commons/PolarizationDIS.h"
      include "../commons/ProjectileDIS.h"
      include "../commons/TargetDIS.h"
      include "../commons/ZedMass.h"
      include "../commons/WMass.h"
      include "../commons/ProtonMass.h"
      include "../commons/SinThetaW.h"
      include "../commons/GFermi.h"
      include "../commons/CKM.h"
*
*     Set all the initialization flags to "xxxx" so that
*     they will be initialized again by InitializeAPFEL.
*
      InWelcome         = "xxxx"
      InScales          = "xxxx"
      InPt              = "xxxx"
      InEvs             = "xxxx"   
      InTheory          = "xxxx"
      InAlpQCD          = "xxxx"
      InAlpQED          = "xxxx"
      InKren            = "xxxx"  
      InMasses          = "xxxx"
      InMFP             = "xxxx"
      InMFA             = "xxxx"
      InPDFs            = "xxxx"
      InRep             = "xxxx"   
      InEvolOp          = "xxxx"
      InLock            = "xxxx"
      InGrid            = "xxxx"
      InTimeLike        = "xxxx"
      InSmallx          = "xxxx"
      InAlphaEvol       = "xxxx"
      InLambdaQCD       = "xxxx"
      InPDFEvol         = "xxxx"
      InFastEvol        = "xxxx"
*
      InMassScheme      = "xxxx"
      inProcessDIS      = "xxxx"
      inPolarizationDIS = "xxxx"
      inProjectileDIS   = "xxxx"
      inTargetDIS       = "xxxx"
*
      InMZ              = "xxxx"
      InMW              = "xxxx"
      InMProton         = "xxxx"
      InSinThetaW       = "xxxx"
      InGFermi          = "xxxx"
      InCKM             = "xxxx"
*
      return
      end
