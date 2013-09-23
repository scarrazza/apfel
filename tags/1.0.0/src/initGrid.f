************************************************************************
*
*     initGrid.f:
*
*     This routine constructs a logarithmically spaced x grid according
*     to the parameters given in input.dat and stores it in a common
*     block.
*
************************************************************************
      subroutine initGrid
*
      implicit none
*
      include "../commons/grid.h"
*
*     Variables
*
      integer ix
*
*     Check
*
      if(igrid.le.0.or.igrid.gt.ngrid)then
         write(6,*) "In initGrid.f:"
         write(6,*) "Invalid value of igrid =",igrid
         call exit(-10)
      endif
*
      do ix=0,nint_max
         xg(igrid,ix) = 0d0
      enddo
*
*     Further check
*
      if(nin(igrid)+inter_degree(igrid).gt.nint_max)then
         write(6,*) "In initGrid.f:"
         write(6,*) "Number of grid points too large:"
         write(6,*) "Maximum value allowed =",nint_max
         write(6,*) "Reduce it in input.dat"
         write(6,*) " "
         call exit(-10)
      endif
*
*     Contruct the baseline grid
*
      step(igrid) = ( log(xmax) - log(xmin(igrid)) ) / dble(nin(igrid))
*
      xg(igrid,0) = xmin(igrid)
      do ix=1,nin(igrid)+inter_degree(igrid)
         xg(igrid,ix) = xg(igrid,ix-1) * exp( step(igrid) )
      enddo
*
      return
      end
