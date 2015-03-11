************************************************************************
*
*     initIntegralsQCD.f:
*
*     This routine initializes the integrals of splitting functions and
*     interpolation functions for a given number of active flavours nf 
*     in QCD.
*
************************************************************************
      subroutine initIntegralsQCD(nf)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/TimeLike.h"
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
*
*     Initialize integrals fo the time-like evolution
*
         if(IsExt(igrid))then
*
*     If this is an external grid, compute the integrals for
*     the entire splitting matrix ...
*
            do alpha=0,nin(igrid)-1
               do beta=alpha,nin(igrid)-1
                  call RSLintegralsQCDT(nf,alpha,beta)
               enddo
            enddo
         else
*
*     ... otherwise only for the first line
*
            do alpha=0,nin(igrid)-1
               call RSLintegralsQCDT(nf,0,alpha)
            enddo
         endif
      else
*
*     Initialize integrals fo the space-like evolution
*
         if(IsExt(igrid))then
*
*     If this is an external grid, compute the integrals for
*     the entire splitting matrix ...
*
            do alpha=0,nin(igrid)-1
               do beta=alpha,nin(igrid)-1
                  call RSLintegralsQCD(nf,alpha,beta)
               enddo
            enddo
         else
*
*     ... otherwise only for the first line
*
            do alpha=0,nin(igrid)-1
               call RSLintegralsQCD(nf,0,alpha)
            enddo
         endif
      endif
*
      return
      end
