************************************************************************
*
*     FLtop.f:
*
*     This function returns the value of the inclusive structure function
*     FL.
*
************************************************************************
      function FLtop(x)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/StructureFunctions.h"
      include "../commons/TMC.h"
      include "../commons/TimeLike.h"
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
      double precision tol
      parameter(tol=1d-10)
**
*     Output Variables
*
      double precision FLtop
*
      if(TMC)then
         tau = 1d0 + 4d0 * rhop * x**2
         xi  = 2d0 * x / ( 1d0 + dsqrt(tau) )
*
         c1 = ( 1d0 - tau ) * x**2 / xi**2 / tau**1.5d0
         c2 = ( 6d0 - 2d0 * tau ) * rhop * x**3 / tau**2
*
         if(xi.lt.xmin(1).or.xi.gt.xmax+tol)then
            write(6,*) "In FLtop.f:"
            write(6,*) "Invalid value of x =",xi
            call exit(-10)
         endif
         if (xi.gt.xmax) xi = 1d0
*
*     Interpolation
*
         FLtop = 0d0
         n = inter_degree(0)
         do alpha=0,nin(0)
            FLtop = FLtop + w_int_gen(n,alpha,xi) 
     1           * ( FL(6,0,alpha) 
     2           + c1 * F2(6,0,alpha) + c2 * I2(6,0,alpha) )
         enddo
         if(dabs(FLtop).le.1d-14) FLtop = 0d0
      else
         if(x.lt.xmin(1).or.x.gt.xmax+tol)then
            write(6,*) "In FLtop.f:"
            write(6,*) "Invalid value of x =",x
            call exit(-10)
         endif
         if (x.gt.xmax) x = 1d0
*
*     Interpolation
*
         FLtop = 0d0
         n = inter_degree(0)
         do alpha=0,nin(0)
            FLtop = FLtop + w_int_gen(n,alpha,x) * FL(6,0,alpha)
         enddo
         if(dabs(FLtop).le.1d-14) FLtop = 0d0
      endif
*
      if(Timelike) FLtop = FLtop / x
*
      return
      end
