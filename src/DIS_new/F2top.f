************************************************************************
*
*     F2top.f:
*
*     This function returns the value of the inclusive structure function
*     F2.
*
************************************************************************
      function F2top(x)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/StructureFunctions.h"
      include "../commons/ProcessDIS.h"
      include "../commons/MassScheme.h"
      include "../commons/m2th.h"
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
      double precision F2top
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In F2top.f:"
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
 101  F2top = 0d0
      n = inter_degree(igrid)
      do alpha=0,nin(igrid)
         F2top = F2top + w_int(n,alpha,x) * F2(6,igrid,alpha)
      enddo
      if(dabs(F2top).le.1d-14) F2top = 0d0
*
      return
      end
