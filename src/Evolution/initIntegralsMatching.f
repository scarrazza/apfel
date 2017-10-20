************************************************************************
*
*     initIntegralsMatching.f:
*
*     This routine initializes the integrals of matching conditions and
*     and interpolation functions.
*
************************************************************************
      subroutine initIntegralsMatching(nf)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/TimeLike.h"
      include "../commons/Smallx.h"
**
*     Input Variables
*
      integer nf
**
*     Internal Variables
*
      integer alpha,beta
*
      if(TimeLike)then
         if(IsExt(igrid))then
            do alpha=0,nin(igrid)-1
               do beta=alpha,nin(igrid)-1
                  call RSLintegralsMatchingT(nf,alpha,beta)
               enddo
            enddo
         else
            do alpha=0,nin(igrid)-1
               call RSLintegralsMatchingT(nf,0,alpha)
            enddo
         endif
      else
         if(IsExt(igrid))then
            do alpha=0,nin(igrid)-1
               do beta=alpha,nin(igrid)-1
                  call RSLintegralsMatching(nf,alpha,beta)
               enddo
            enddo
            if(Smallx.and.nf.lt.6)then
               do alpha=0,nin(igrid)-1
                  do beta=alpha,nin(igrid)-1
                     call RSLintegralsMatchingRes(nf,alpha,beta)
                  enddo
               enddo
            endif
         else
            do alpha=0,nin(igrid)-1
               call RSLintegralsMatching(nf,0,alpha)
            enddo
            if(Smallx.and.nf.lt.6)then
               do alpha=0,nin(igrid)-1
                  call RSLintegralsMatchingRes(nf,0,alpha)
               enddo
            endif
         endif
      endif
*
      return
      end
