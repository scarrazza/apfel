************************************************************************
*
*     initIntegralsQED.f:
*
*     This routine initializes the integrals of splitting functions and
*     and interpolation functions for a given number of active flavours 
*     nf in QCD.
*
************************************************************************
      subroutine initIntegralsQED(nf)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      integer nf
**
*     Internal Variables
*
      integer alpha
*
      do alpha=0,nin(igrid)-1
         call RSLintegralsQED(nf,0,alpha)
      enddo
*
      return
      end
