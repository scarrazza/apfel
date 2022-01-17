************************************************************************
*
*     g1charm.f:
*
*     This function returns the value of the charm structure function
*     g1c.
*
************************************************************************
      function g1charm(x)
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
**
*     Output Variables
*
      double precision g1charm
*
      if(Polarized)then
*      
         if(TMC)then
            write(6,*) "TMCs not available for polarised DIS"
            call exit(-10)      
         else
            if(x.lt.xmin(1)-tol.or.x.gt.xmax+tol)then
               write(6,*) "In g1charm.f:"
               write(6,*) "Invalid value of x =",x
               call exit(-10)
            endif
            if (x.lt.xmin(1)) x = xmin(1)
            if (x.gt.xmax) x = 1d0
*     
*     Interpolation
*     
            g1charm = 0d0
            n = inter_degree(0)
            do alpha=0,nin(0)
               g1charm = g1charm + w_int_gen(n,alpha,x) * F2(4,0,alpha)
            enddo
            if(dabs(g1charm).le.1d-14) g1charm = 0d0
         endif
      else
         write(6,*) "g1 structure function not available",
     1        " for unpolarised DIS"
         call exit(-10)
      endif
*     
      return
      end
      
