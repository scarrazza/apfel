************************************************************************
*
*     F2total.f:
*
*     This function returns the value of the inclusive structure function
*     F2.
*
************************************************************************
      function F2total(x)
*
      implicit none
**
*     Input Variables
*
      double precision x
**
*     Internal Variables
*
      double precision F2light,F2charm,F2bottom,F2top
**
*     Output Variables
*
      double precision F2total
*
      F2total = F2light(x)
     1        + F2charm(x)
     2        + F2bottom(x)
     3        + F2top(x)
*
      return
      end
c$$$*
c$$$************************************************************************
c$$$      function F2total(x)
c$$$*
c$$$      implicit none
c$$$*
c$$$      include "../commons/grid.h"
c$$$      include "../commons/StructureFunctions.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      double precision x
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      integer n
c$$$      integer alpha
c$$$      double precision w_int
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      double precision F2total
c$$$*
c$$$      if(x.lt.xmin(1).or.x.gt.xmax)then
c$$$         write(6,*) "In F2total.f:"
c$$$         write(6,*) "Invalid value of x =",x
c$$$         call exit(-10)
c$$$      endif
c$$$*
c$$$*     Select the grid
c$$$*
c$$$      do igrid=1,ngrid
c$$$         if(x.ge.xmin(igrid).and.x.lt.xmin(igrid+1))then
c$$$            goto 101
c$$$         endif
c$$$      enddo
c$$$*
c$$$*     Interpolation
c$$$*
c$$$ 101  F2total = 0d0
c$$$      n = inter_degree(igrid)
c$$$      do alpha=0,nin(igrid)
c$$$         F2total = F2total + w_int(n,alpha,x) * F2(7,igrid,alpha)
c$$$      enddo
c$$$      if(dabs(F2total).le.1d-14) F2total = 0d0
c$$$*
c$$$      return
c$$$      end
