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
      integer alpha,beta
*
      if(TimeLike)then
         if(IsExt(igrid))then
            do alpha=0,nin(igrid)-1
               do beta=alpha,nin(igrid)-1
                  call RSLintegralsMatchingT(alpha,beta)
               enddo
            enddo
         else
            do alpha=0,nin(igrid)-1
               call RSLintegralsMatchingT(0,alpha)
            enddo
         endif
      else
         if(IsExt(igrid))then
            do alpha=0,nin(igrid)-1
               do beta=alpha,nin(igrid)-1
                  call RSLintegralsMatching(alpha,beta)
               enddo
            enddo
         else
            do alpha=0,nin(igrid)-1
               call RSLintegralsMatching(0,alpha)
            enddo
         endif
      endif
*
      return
      end
