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
      if(inMProton.ne."done")then
         write(6,*) "GetProtonMass: Parameter not initialized"
         write(6,*) "Set it by means of 'SetProtonMass'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      GetProtonMass = MProton
*
      return
      end
