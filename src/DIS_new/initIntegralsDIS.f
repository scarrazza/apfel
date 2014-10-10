************************************************************************
*
*     initIntegralsDIS.f:
*
*     This routine initializes the integrals of coefficient functions
*     and interpolation functions.
*
************************************************************************
      subroutine initIntegralsDIS
*
      implicit none
*
      include "../commons/grid.h"
**
*     Internal Variables
*
      integer alpha,beta
*
*     Initialize integrals
*
      if(IsExt(igrid))then
*
*     If this is an external grid, compute the integrals for
*     the entire splitting matrix ...
*
         do alpha=0,nin(igrid)-1
            do beta=alpha,nin(igrid)-1
               call RSLintegralsDIS(alpha,beta)
            enddo
         enddo
      else
*
*     ... otherwise only for the first line
*
         do alpha=0,nin(igrid)-1
            call RSLintegralsDIS(0,alpha)
         enddo
      endif
*
      return
      end
