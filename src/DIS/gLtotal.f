************************************************************************
*
*     gLtotal.f:
*
*     This function returns the value of the inclusive structure function
*     gL.
*
************************************************************************
      function gLtotal(x)
*
      implicit none
**
*     Input Variables
*
      double precision x
**
*     Internal Variables
*
      double precision gLlight,gLcharm,gLbottom,gLtop
**
*     Output Variables
*
      double precision gLtotal
*
      gLtotal = gLlight(x)
     1        + gLcharm(x)
     2        + gLbottom(x)
     3        + gLtop(x)
*
      return
      end
*
c$$$************************************************************************
c$$$      function FLtotal(x)
c$$$*
c$$$      implicit none
c$$$*
c$$$      include "../commons/grid.h"
c$$$      include "../commons/StructureFunctions.h"
c$$$      include "../commons/TMC.h"
c$$$      include "../commons/TimeLike.h"
c$$$**
c$$$*     Input Variables
c$$$*
c$$$      double precision x
c$$$**
c$$$*     Internal Variables
c$$$*
c$$$      integer n
c$$$      integer alpha
c$$$      double precision w_int_gen
c$$$      double precision tau,xi
c$$$      double precision c1,c2
c$$$      double precision tol
c$$$      parameter(tol=1d-10)
c$$$**
c$$$*     Output Variables
c$$$*
c$$$      double precision FLtotal
c$$$*
c$$$      if(TMC)then
c$$$         tau = 1d0 + 4d0 * rhop * x**2
c$$$         xi  = 2d0 * x / ( 1d0 + dsqrt(tau) )
c$$$*
c$$$         c1 = ( 1d0 - tau ) * x**2 / xi**2 / tau**1.5d0
c$$$         c2 = ( 6d0 - 2d0 * tau ) * rhop * x**3 / tau**2
c$$$*
c$$$         if(xi.lt.xmin(1)-tol.or.xi.gt.xmax+tol)then
c$$$            write(6,*) "In FLtotal.f:"
c$$$            write(6,*) "Invalid value of x =",xi
c$$$            call exit(-10)
c$$$         endif
c$$$         if (xi.lt.xmin(1)) xi = xmin(1)
c$$$         if (xi.gt.xmax) xi = 1d0
c$$$*
c$$$*     Interpolation
c$$$*
c$$$         FLtotal = 0d0
c$$$         n = inter_degree(0)
c$$$         do alpha=0,nin(0)
c$$$            FLtotal = FLtotal + w_int_gen(n,alpha,xi) 
c$$$     1           * ( FL(7,0,alpha) 
c$$$     2           + c1 * F2(7,0,alpha) + c2 * I2(7,0,alpha) )
c$$$         enddo
c$$$         if(dabs(FLtotal).le.1d-14) FLtotal = 0d0
c$$$      else
c$$$         if(x.lt.xmin(1)-tol.or.x.gt.xmax+tol)then
c$$$            write(6,*) "In FLtotal.f:"
c$$$            write(6,*) "Invalid value of x =",x
c$$$            call exit(-10)
c$$$         endif
c$$$         if (x.lt.xmin(1)) x = xmin(1)
c$$$         if (x.gt.xmax) x = 1d0
c$$$*
c$$$*     Interpolation
c$$$*
c$$$         FLtotal = 0d0
c$$$         n = inter_degree(0)
c$$$         do alpha=0,nin(0)
c$$$            FLtotal = FLtotal + w_int_gen(n,alpha,x) * FL(7,0,alpha)
c$$$         enddo
c$$$         if(dabs(FLtotal).le.1d-14) FLtotal = 0d0
c$$$      endif
c$$$*
c$$$      if(Timelike) FLtotal = FLtotal / x
c$$$*
c$$$      return
c$$$      end
