************************************************************************
*
*     g4total.f:
*
*     This function returns the value of the inclusive structure function
*     g4.
*
************************************************************************
      function g4total(x)
*
      implicit none
**
*     Input Variables
*
      double precision x
**
*     Internal Variables
*
      double precision g4light,g4charm,g4bottom,g4top
**
*     Output Variables
*
      double precision g4total
*
      g4total = g4light(x)
     1        + g4charm(x)
     2        + g4bottom(x)
     3        + g4top(x)
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
c$$$      double precision F2total
c$$$*
c$$$      if(TMC)then
c$$$         tau = 1d0 + 4d0 * rhop * x**2
c$$$         xi  = 2d0 * x / ( 1d0 + dsqrt(tau) )
c$$$*
c$$$         c1 = x**2 / xi**2 / tau**1.5d0
c$$$         c2 = 6d0 * rhop * x**3 / tau**2
c$$$*
c$$$         if(xi.lt.xmin(1)-tol.or.xi.gt.xmax+tol)then
c$$$            write(6,*) "In F2total.f:"
c$$$            write(6,*) "Invalid value of x =",xi
c$$$            call exit(-10)
c$$$         endif
c$$$         if (xi.lt.xmin(1)) xi = xmin(1)
c$$$         if (xi.gt.xmax) xi = 1d0
c$$$*
c$$$*     Interpolation
c$$$*
c$$$         F2total = 0d0
c$$$         n = inter_degree(0)
c$$$         do alpha=0,nin(0)
c$$$            F2total = F2total + w_int_gen(n,alpha,xi) 
c$$$     1           * ( c1 * F2(7,0,alpha) + c2 * I2(7,0,alpha) )
c$$$         enddo
c$$$         if(dabs(F2total).le.1d-14) F2total = 0d0
c$$$      else
c$$$         if(x.lt.xmin(1)-tol.or.x.gt.xmax+tol)then
c$$$            write(6,*) "In F2total.f:"
c$$$            write(6,*) "Invalid value of x =",x
c$$$            call exit(-10)
c$$$         endif
c$$$         if (x.lt.xmin(1)) x = xmin(1)
c$$$         if (x.gt.xmax) x = 1d0
c$$$*
c$$$*     Interpolation
c$$$*
c$$$         F2total = 0d0
c$$$         n = inter_degree(0)
c$$$         do alpha=0,nin(0)
c$$$            F2total = F2total + w_int_gen(n,alpha,x) * F2(7,0,alpha)
c$$$         enddo
c$$$         if(dabs(F2total).le.1d-14) F2total = 0d0
c$$$      endif
c$$$*
c$$$      if(Timelike) F2total = F2total / x
c$$$*
c$$$      return
c$$$      end
