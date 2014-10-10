************************************************************************
*
*     FLlight.f:
*
*     This function returns the value of the light structure function
*     FLl.
*
************************************************************************
      function FLlight(x)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/StructureFunctions.h"
**
*     Input Variables
*
      double precision x
**
*     Internal Variables
*
      integer n
      integer alpha
      integer i
      double precision w_int
**
*     Output Variables
*
      double precision FLlight
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In FLlight.f:"
         write(6,*) "Invalid value of x =",x
         call exit(-10)
      endif
*
*     Select the grid
*
      do igrid=1,ngrid
         if(x.ge.xmin(igrid).and.x.lt.xmin(igrid+1))then
            goto 101
         endif
      enddo
*
*     Interpolation
*
 101  FLlight = 0d0
      n = inter_degree(igrid)
      do alpha=0,nin(igrid)
         do i=1,3
            FLlight = FLlight + w_int(n,alpha,x) * FL(i,igrid,alpha)
         enddo
      enddo
      if(dabs(FLlight).le.1d-14) FLlight = 0d0
*
      return
      end
