************************************************************************
*
*     F3total.f:
*
*     This function returns the value of the inclusive structure function
*     F3.
*
************************************************************************
      function F3total(x)
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
      double precision w_int
**
*     Output Variables
*
      double precision F3total
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In F3total.f:"
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
 101  F3total = 0d0
      n = inter_degree(igrid)
      do alpha=0,nin(igrid)
         F3total = F3total + w_int(n,alpha,x) * F3(7,igrid,alpha)
      enddo
      if(dabs(F3total).le.1d-14) F3total = 0d0
*
      return
      end
