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
      if(inMW.ne."done")then
         write(6,*) "GetWMass: Parameter not initialized"
         write(6,*) "Set it by means of 'SetWMass'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      GetWMass = MW
*
      return
      end
