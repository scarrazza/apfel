************************************************************************
*
*     GetMaxFlavourAlpha.f:
*
*     This function gets the maximum number of flavours that the
*     evolution of alphaQCD and alphaQED can reach.
*
************************************************************************
      function GetMaxFlavourAlpha()
*
      implicit none
*
      include "../commons/MaxFlavourAlpha.h"
*
*     Variables
*
      integer GetMaxFlavourAlpha
*
      if(InMFA.ne."done")then
         write(6,*) "GetMaxFlavourAlpha: Parameter not initialized"
         write(6,*) "Set it by means of 'SetMaxFlavourAlpha'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      GetMaxFlavourAlpha = nfMaxAlpha
*
      return
      end
