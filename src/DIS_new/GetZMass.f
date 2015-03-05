************************************************************************
*
*     GetZMass.f:
*
*     This function returns the mass of the Z in GeV.
*
************************************************************************
      function GetZMass()
*
      implicit none
*
      include "../commons/ZedMass.h"
*
*     Variables
*
      double precision GetZMass
*
      GetZMass = MZ
*
      return
      end
