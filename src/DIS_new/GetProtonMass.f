************************************************************************
*
*     GetProtonMass.f:
*
*     This function returns the mass of the Proton in GeV.
*
************************************************************************
      function GetProtonMass()
*
      implicit none
*
      include "../commons/ProtonMass.h"
*
*     Variables
*
      double precision GetProtonMass
*
      GetProtonMass = MProton
*
      return
      end
