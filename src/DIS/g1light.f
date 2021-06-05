************************************************************************
*
*     g1light.f:
*
*     This function returns the value of the light structure function
*     g1l.
*
************************************************************************
      function g1light(x)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/StructureFunctions.h"
      include "../commons/TMC.h"
      include "../commons/Polarized.h"
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
      double precision tol
      parameter(tol=1d-10)
      double precision F3light
**
*     Output Variables
*
      double precision g1light
*
      if(Polarized)then
*      
         if(TMC)then
            write(6,*) "TMCs not available for polarised DIS"
            call exit(-10)
         else
            if(x.lt.xmin(1)-tol.or.x.gt.xmax+tol)then
               write(6,*) "In g1light.f:"
               write(6,*) "Invalid value of x =",x
               call exit(-10)
            endif
            if (x.lt.xmin(1)) x = xmin(1)
            if (x.gt.xmax) x = 1d0
*     
*     Interpolation
*     
            F3light = 0d0
            n = inter_degree(0)
            do alpha=0,nin(0)
               F3light = F3light + w_int_gen(n,alpha,x) * F3(3,0,alpha)
            enddo
            if(dabs(F3light).le.1d-14) F3light = 0d0
         endif
*
         g1light = F3light / 2d0
*
      else
         write(6,*) "g1 structure function not available",
     1        " for unpolarised DIS"
         call exit(-10)
      endif
* 
      return
      end
