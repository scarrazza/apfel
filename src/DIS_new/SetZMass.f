************************************************************************
*
*     SetZMass.f:
*
*     This subroutine sets the mass of the Z in GeV.
*
************************************************************************
      subroutine SetZMass(massz)
*
      implicit none
*
      include "../commons/ZedMass.h"
*
*     Variables
*
      double precision massz
*
      MZ   = massz
      InMZ = "done"
*
      return
      end
