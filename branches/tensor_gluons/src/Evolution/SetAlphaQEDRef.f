************************************************************************
*
*     SetAlphaQEDRef.f:
*
*     This subroutine sets the reference value of alpha at the
*     reference scale.
*
************************************************************************
      subroutine SetAlphaQEDRef(alpharef,Qref)
*
      implicit none
*
      include "../commons/alpha_ref_QED.h"
*
*     Variables
*
      double precision alpharef,Qref
*
      alpha_ref_QED = alpharef
      q2_ref_QED    = Qref * Qref
      InAlpQED      = "done"
*
      return
      end
