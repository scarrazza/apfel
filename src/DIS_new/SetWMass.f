************************************************************************
*
*     SetWMass.f:
*
*     This subroutine sets the mass of the W in GeV.
*
************************************************************************
      subroutine SetWMass(massw)
*
      implicit none
*
      include "../commons/WMass.h"
*
*     Variables
*
      double precision massw
*
      MW   = massw
      InMW = "done"
*
      return
      end
