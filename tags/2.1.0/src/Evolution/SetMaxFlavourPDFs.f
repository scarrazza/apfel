************************************************************************
*
*     SetMaxFlavourPDFs.f:
*
*     This subroutine sets the maximum number of flavours that the
*     evolution of PDFs can reach.
*
************************************************************************
      subroutine SetMaxFlavourPDFs(nf)
*
      implicit none
*
      include "../commons/MaxFlavourPDFs.h"
*
*     Variables
*
      integer nf
*
      nfMaxPDFs = nf
      InMFP     = "done"
*
      return
      end
