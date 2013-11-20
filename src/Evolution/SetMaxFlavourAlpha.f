************************************************************************
*
*     SetMaxFlavourAlpha.f:
*
*     This subroutine sets the maximum number of flavours that the
*     evolution of alphaQCD and alphaQED can reach.
*
************************************************************************
      subroutine SetMaxFlavourAlpha(nf)
*
      implicit none
*
      include "../commons/MaxFlavourAlpha.h"
*
*     Variables
*
      integer nf
*
      nfMaxAlpha = nf
      InMFA      = "done"
*
      return
      end
