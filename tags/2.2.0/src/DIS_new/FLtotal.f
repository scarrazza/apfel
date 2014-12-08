************************************************************************
*
*     FLtotal.f:
*
*     This function returns the value of the inclusive structure function
*     FL.
*
************************************************************************
      function FLtotal(x)
*
      implicit none
**
*     Input Variables
*
      double precision x
**
*     Internal Variables
*
      double precision FLlight,FLcharm,FLbottom,FLtop
**
*     Output Variables
*
      double precision FLtotal
*
      FLtotal = FLlight(x)
     1        + FLcharm(x)
     2        + FLbottom(x)
     3        + FLtop(x)
*
      return
      end
c$$$*
c$$$************************************************************************
c$$$      function FLtotal(x)
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
c$$$      double precision FLtotal
c$$$*
c$$$      if(x.lt.xmin(1).or.x.gt.xmax)then
c$$$         write(6,*) "In FLtotal.f:"
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
c$$$ 101  FLtotal = 0d0
c$$$      n = inter_degree(igrid)
c$$$      do alpha=0,nin(igrid)
c$$$         FLtotal = FLtotal + w_int(n,alpha,x) * FL(7,igrid,alpha)
c$$$      enddo
c$$$      if(dabs(FLtotal).le.1d-14) FLtotal = 0d0
c$$$*
c$$$      return
c$$$      end
