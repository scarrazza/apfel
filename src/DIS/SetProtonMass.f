************************************************************************
*
*     SetProtonMass.f:
*
*     This subroutine sets the mass of the Proton in GeV.
*
************************************************************************
      subroutine SetProtonMass(massp)
*
      implicit none
*
      include "../commons/ProtonMass.h"
*
*     Variables
*
      double precision massp
*
      MProton   = massp
      InMProton = "done"
*
      return
      end
