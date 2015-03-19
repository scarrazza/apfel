************************************************************************
*
*     SetTauMass.f:
*
*     This subroutine sets the mass of the Tau in GeV.
*
************************************************************************
      subroutine SetTauMass(masst)
*
      implicit none
*
      include "../commons/TauMass.h"
*
*     Variables
*
      double precision masst
*
      MTau   = masst
      InMTau = "done"
*
      return
      end
