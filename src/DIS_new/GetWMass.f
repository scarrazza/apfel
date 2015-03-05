************************************************************************
*
*     GetWMass.f:
*
*     This function returns the mass of the W in GeV.
*
************************************************************************
      function GetWMass()
*
      implicit none
*
      include "../commons/WMass.h"
*
*     Variables
*
      double precision GetWMass
*
      GetWMass = MW
*
      return
      end
