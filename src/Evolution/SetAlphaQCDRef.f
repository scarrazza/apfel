************************************************************************
*
*     SetAlphaQCDRef.f:
*
*     This subroutine sets the reference value of alpha_s at the
*     reference scale.
*
************************************************************************
      subroutine SetAlphaQCDRef(alpharef,Qref)
*
      implicit none
*
      include "../commons/alpha_ref_QCD.h"
      include "../commons/InAPFEL.h"
*
*     Variables
*
      double precision alpharef,Qref
*
      alpha_ref_QCD = alpharef
      q2_ref_QCD    = Qref * Qref
      InAlpQCD      = "done"
*
*     Upadate the values of alphas at the thresholds.
*     Do it only if this function is called after APFEL
*     has been initialized.
*
      if(InAPFEL.eq."done") call ThresholdAlphaQCD
*
      return
      end
