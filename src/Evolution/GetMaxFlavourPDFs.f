************************************************************************
*
*     GetMaxFlavourPDFs.f:
*
*     This function gets the maximum number of flavours that the
*     evolution of PDFs can reach.
*
************************************************************************
      function GetMaxFlavourPDFs()
*
      implicit none
*
      include "../commons/MaxFlavourPDFs.h"
*
*     Variables
*
      integer GetMaxFlavourPDFs
*
      if(InMFP.ne."done")then
         write(6,*) "GetMaxFlavourPDFs: Parameter not initialized"
         write(6,*) "Set it by means of 'SetMaxFlavourPDFs'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      GetMaxFlavourPDFs = nfMaxPDFs
*
      return
      end
