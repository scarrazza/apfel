************************************************************************
*
*     initIntegralsQED.f:
*
*     This routine initializes the integrals of splitting functions and
*     interpolation functions for a given number of active flavours nf 
*     in QCD.
*
************************************************************************
      subroutine initIntegralsQED(nf,nl)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      integer nf,nl
**
*     Internal Variables
*
      integer alpha
*
      do alpha=0,nin(igrid)-1
         call RSLintegralsQED(nf,nl,0,alpha)
      enddo
*
      return
      end
