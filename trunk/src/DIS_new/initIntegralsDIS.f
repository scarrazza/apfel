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
*     Check that the number of intervals does not exceed the maximum number
*
      if(nin(igrid)+inter_degree(igrid).gt.nint_max_DIS)then
         write(6,*) "In initIntegralsDIS.f:"
         write(6,*) "Number of grid points too large:"
         write(6,*) "Maximum value allowed =",nint_max_DIS
         write(6,*) "You should reduce it."
         write(6,*) " "
         call exit(-10)
      endif
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
