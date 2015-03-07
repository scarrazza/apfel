************************************************************************
*
*     GetGFermi.f:
*
*     This function returns the Fermi constant.
*
************************************************************************
      function GetGFermi()
*
      implicit none
*
      include "../commons/GFermi.h"
*
*     Variables
*
      double precision GetGFermi
*
      if(inGFermi.ne."done")then
         write(6,*) "GetGFermi: Parameter not initialized"
         write(6,*) "Set it by means of 'SetGFermi'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      GetGFermi = GFermi
*
      return
      end
