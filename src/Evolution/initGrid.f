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
      integer alpha
      integer jgrid
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
*     Now construct the grid
*
      if(IsExt(igrid))then
         do ix=0,nin(igrid)
            xg(igrid,ix) = xgext(igrid,ix)
         enddo
*
*     Check that the last point of the user-given grid is equal to one
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
         xg(igrid,nin(igrid)) = 1d0
      endif
*
*     Join grids (do it only after the last grid has been computed)
*
      if(igrid.eq.ngrid)then
         nin(0) = -1
         TransitionPoint(1) = 0
         do jgrid=1,ngrid
            do alpha=0,nin(jgrid)
               if(xmin(jgrid+1)-xg(jgrid,alpha).lt.eps)then
                  TransitionPoint(jgrid+1) = nin(0) + 1
                  goto 101
               endif
               nin(0) = nin(0) + 1
*     Joining x-space grid ...
               xg(0,nin(0)) = xg(jgrid,alpha)
            enddo
 101     enddo
         TransitionPoint(ngrid+1) = nin(0)
*
*     Ensure that the last point of the joint grid is
*     exactly equal to 1.
*
         if(dabs(xg(0,nin(0))-1d0).gt.eps) xg(0,nin(0)) = 1d0
*
*     (Indicative) interpolation degree of the joint grid
*     just to be used in the check below.
*
         inter_degree(0) = inter_degree(ngrid)
c         inter_degree(0) = inter_degree(1)
*
*     Grid for x > 1 (Needed by the interpolation)
*
         do alpha=1,inter_degree(0)
            xg(0,nin(0)+alpha) = xg(ngrid,nin(ngrid)+alpha)
         enddo
*
         if(nin(0)+inter_degree(0).gt.nint_max)then
            write(6,*) "In initGrids.f:"
            write(6,*) "Number of points of the joint grid too large:"
            write(6,*) "Maximum value allowed =",nint_max
            write(6,*) "found =",nin(0)+inter_degree(0)
            write(6,*) "You should reduce it"
            write(6,*) " "
            call exit(-10)
         endif
      endif
*
      return
      end
