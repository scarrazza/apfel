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
      include "../commons/TimeLike.h"
**
*     Internal Variables
*
      integer alpha
*
      if(TimeLike)then
         do alpha=0,nin(igrid)-1
            call RSLintegralsMatchingT(0,alpha)
         enddo
      else
         do alpha=0,nin(igrid)-1
            call RSLintegralsMatching(0,alpha)
         enddo
      endif
*
      return
      end
