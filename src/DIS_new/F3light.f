************************************************************************
*
*     F3light.f:
*
*     This function returns the value of the light structure function
*     F3l.
*
************************************************************************
      function F3light(x)
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
      double precision F3light
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In F3light.f:"
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
 101  F3light = 0d0
      n = inter_degree(igrid)
      do alpha=0,nin(igrid)
         do i=1,3
            F3light = F3light + w_int(n,alpha,x) * F3(i,igrid,alpha)
         enddo
      enddo
      if(dabs(F3light).le.1d-14) F3light = 0d0
*
      return
      end