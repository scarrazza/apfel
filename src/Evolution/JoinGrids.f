************************************************************************
*
*     JoinGrids.f:
*
************************************************************************
      subroutine JoinGrids
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/fph.h"
**
*     Internal Variables
*
      integer alpha,jgrid,i
      double precision eps
      parameter(eps=1d-20)
*
*     Join grids
*
      nin(0) = -1
      do jgrid=1,ngrid
         do alpha=0,nin(jgrid)
            if(xmin(jgrid+1)-xg(jgrid,alpha).lt.eps)then
               goto 102
            endif
            nin(0) = nin(0) + 1
*     Joining x-space grid ...
            xg(0,nin(0)) = xg(jgrid,alpha)
*     Joining PDFs ...
            do i=-6,6
               fph(0,i,nin(0)) = fph(jgrid,i,alpha)
            enddo
            fgamma(0,nin(0)) = fgamma(jgrid,alpha)
         enddo
 102  enddo
*
*     (Indicative) interpolation degree of the joint grid
*     just to be used in the check below.
*
      inter_degree(0) = inter_degree(ngrid)
*
*     Grid and PDFs for x > 1 (Needed by the interpolation)
*
      do alpha=1,inter_degree(0)
         xg(0,nin(0)+alpha) = xg(ngrid,nin(ngrid)+alpha)
         do i=-6,6
            fph(0,i,nin(0)+alpha) = fph(ngrid,i,nin(ngrid)+alpha)
         enddo
         fgamma(0,nin(0)+alpha) = fgamma(ngrid,nin(ngrid)+alpha)
      enddo
*
      if(nin(0)+inter_degree(0).gt.nint_max)then
         write(6,*) "In JoinGrids.f:"
         write(6,*) "Number of points of the joint grid too large:"
         write(6,*) "Maximum value allowed =",nint_max
         write(6,*) "You should reduce it"
         write(6,*) " "
         call exit(-10)
      endif
*
      open(unit=19,file="JointGrid.dat",status="unknown")
      write(19,*) nin(0)!+inter_degree(0)
      do alpha=0,nin(0)!+inter_degree(0)
         write(19,*) xg(0,alpha)
      enddo
      close(19)
*
      return
      end
