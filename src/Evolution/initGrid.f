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
      include "../commons/lock.h"
      include "../commons/grid.h"
*     
*     Variables
*     
      integer ix
      double precision eps
      parameter(eps=1d-12)
*
*     Checks
*
      if(igrid.le.0.or.igrid.gt.ngrid)then
         write(6,*) "In initGrid.f:"
         write(6,*) "Invalid value of igrid =",igrid
         call exit(-10)
      endif
*
      if(nin(igrid)+inter_degree(igrid).gt.nint_max)then
         write(6,*) "In initGrid.f:"
         write(6,*) "Number of grid points too large:"
         write(6,*) "found =",nin(igrid)+inter_degree(igrid)
         write(6,*) "Maximum value allowed =",nint_max
         write(6,*) "You should reduce it."
         write(6,*) " "
         call exit(-10)
      endif
*
*     Initialize subgrid
*
      do ix=0,nint_max
         xg(igrid,ix) = 0d0
      enddo
*
*     Set density factor to one by default
*
      DensityFactor(igrid) = 1
*
*     In case the grids have been locked ...
*     (Needed to export the evolution operator)
*
      if(lock.and.igrid.gt.1)then
*
*     Find the closest point of the "(igrid-1)"-th subgris to "xmin(igrid)"
*     and replace "xmin(igrid)".
*
         do ix=0,nin(igrid-1)
            if((xmin(igrid)-xg(igrid-1,ix)).le.0d0) goto 102
         enddo
 102     if(dabs(xmin(igrid)-xg(igrid-1,ix)).gt.
     1      dabs(xmin(igrid)-xg(igrid-1,ix+1)))then
            xmin(igrid) = xg(igrid-1,ix+1)
         else
            xmin(igrid) = xg(igrid-1,ix)
         endif
*
*     Find the closest multiple of "nin(igrid-1) - ix + 1" to "nin(igrid)",
*     i.e. "DensityFactor(igrid)", and replace "nin(igrid)"
*
         DensityFactor(igrid) = nint( dble( nin(igrid) ) 
     1                        / dble( nin(igrid-1) - ix ) )
         nin(igrid) = DensityFactor(igrid) * ( nin(igrid-1) - ix )
      endif
*
*     Now contruct the grid
*
      if(IsExt(igrid))then
         do ix=0,nin(igrid)
            xg(igrid,ix) = xgext(igrid,ix)
         enddo
*
*     Check that the last point of the user give grid is is euqal to one
*
         if(dabs(xg(igrid,nin(igrid))-1d0).gt.eps)then
            write(6,*) "In initGrid.f:"
            write(6,*) "The upper bound of the ",igrid,"-th grid",
     1                 " does not coincide with one: xmax =",
     2                 xg(igrid,nin(igrid))
            write(6,*) "Check the input grid"
            call exit(-10)
         else
            xg(igrid,nin(igrid)) = 1d0
         endif
*
*     Extend the grid for x > 1 for interpolation reasons using the same
*     width of the last bin in log scale
*
         step(igrid) = dlog(xg(igrid,nin(igrid))/xg(igrid,nin(igrid)-1))
         do ix=nin(igrid)+1,nin(igrid)+inter_degree(igrid)
            xg(igrid,ix) = xg(igrid,ix-1) * exp( step(igrid) )
         enddo
      else
         step(igrid) = ( dlog(xmax) - dlog(xmin(igrid)) )
     1               / dble(nin(igrid))
*
         xg(igrid,0) = xmin(igrid)
         do ix=1,nin(igrid)+inter_degree(igrid)
            xg(igrid,ix) = xg(igrid,ix-1) * exp( step(igrid) )
         enddo
         xg(igrid, nin(igrid)) = 1d0
      endif
*
      return
      end
