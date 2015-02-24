************************************************************************
*
*     F2light.f:
*
*     This function returns the value of the light structure function
*     F2l.
*
************************************************************************
      function F2light(x)
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
      double precision F2light
*
      if(TMC)then
         tau = 1d0 + 4d0 * rhop * x**2d0
         xi  = 2d0 * x / ( 1d0 + dsqrt(tau) )
*
         c1 = x**2d0 / xi**2d0 / tau**1.5d0
         c2 = 6d0 * rhop * x**3d0 / tau**2d0
*
         if(xi.lt.xmin(1).or.xi.gt.xmax)then
            write(6,*) "In F2light.f:"
            write(6,*) "Invalid value of x =",xi
            call exit(-10)
         endif
*
*     Interpolation
*
         F2light = 0d0
         n = inter_degree(0)
         do alpha=0,nin(0)
            F2light = F2light + w_int_gen(n,alpha,xi) 
     1           * ( c1 * F2(3,0,alpha) + c2 * I2(3,0,alpha) )
         enddo
         if(dabs(F2light).le.1d-14) F2light = 0d0
      else
         if(x.lt.xmin(1).or.x.gt.xmax)then
            write(6,*) "In F2light.f:"
            write(6,*) "Invalid value of x =",x
            call exit(-10)
         endif
*
*     Interpolation
*
         F2light = 0d0
         n = inter_degree(0)
         do alpha=0,nin(0)
            F2light = F2light + w_int_gen(n,alpha,x) * F2(3,0,alpha)
         enddo
         if(dabs(F2light).le.1d-14) F2light = 0d0
      endif
*
      return
      end
