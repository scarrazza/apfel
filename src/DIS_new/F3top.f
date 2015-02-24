************************************************************************
*
*     F3top.f:
*
*     This function returns the value of the inclusive structure function
*     F3.
*
************************************************************************
      function F3top(x)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/StructureFunctions.h"
      include "../commons/TMC.h"
**
*     Input Variables
*
      double precision x
**
*     Internal Variables
*
      integer n
      integer alpha
      double precision w_int_gen
      double precision tau,xi
      double precision c1,c2
**
*     Output Variables
*
      double precision F3top
*
*
      if(TMC)then
         tau = 1d0 + 4d0 * rhop * x**2d0
         xi  = 2d0 * x / ( 1d0 + dsqrt(tau) )
*
         c1 = x**2d0 / xi**2d0 / tau
         c2 = 4d0 * rhop * x**3d0 / tau**1.5d0
*
         if(xi.lt.xmin(1).or.xi.gt.xmax)then
            write(6,*) "In F3top.f:"
            write(6,*) "Invalid value of x =",xi
            call exit(-10)
         endif
*
*     Interpolation
*
         F3top = 0d0
         n = inter_degree(0)
         do alpha=0,nin(0)
            F3top = F3top + w_int_gen(n,alpha,xi) 
     1           * ( c1 * F3(6,0,alpha) + c2 * I3(6,0,alpha) )
         enddo
         if(dabs(F3top).le.1d-14) F3top = 0d0
      else
         if(x.lt.xmin(1).or.x.gt.xmax)then
            write(6,*) "In F3top.f:"
            write(6,*) "Invalid value of x =",x
            call exit(-10)
         endif
*
*     Interpolation
*
         F3top = 0d0
         n = inter_degree(0)
         do alpha=0,nin(0)
            F3top = F3top + w_int_gen(n,alpha,x) * F3(6,0,alpha)
         enddo
         if(dabs(F3top).le.1d-14) F3top = 0d0
      endif
*
      return
      end
