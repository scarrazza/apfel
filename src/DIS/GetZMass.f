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
      if(inMZ.ne."done")then
         write(6,*) "GetZMass: Parameter not initialized"
         write(6,*) "Set it by means of 'SetZMass'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      GetZMass = MZ
*
      return
      end
