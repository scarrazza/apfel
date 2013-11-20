************************************************************************
*
*     SetAlphaQCDRef.f:
*
*     This subroutine sets the reference vales of alpha_s at the
*     reference scale.
*
************************************************************************
      subroutine SetAlphaQCDRef(alpharef,Qref)
*
      implicit none
*
      include "../commons/alpha_ref_QCD.h"
*
*     Variables
*
      double precision alpharef,Qref
*
      alpha_ref_QCD = alpharef
      q2_ref_QCD    = Qref * Qref
      InAlpQCD      = "done"
*
      return
      end
