************************************************************************
*
*     initIntegralsMatching.f:
*
*     This routine initializes the integrals of matching conditions and
*     and interpolation functions.
*
************************************************************************
      subroutine initIntegralsMatching
*
      implicit none
*
      include "../commons/grid.h"
**
*     Internal Variables
*
      integer alpha
*
      do alpha=0,nin(igrid)-1
         call RSLintegralsMatching(0,alpha)
      enddo
*
      return
      end
