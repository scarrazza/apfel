************************************************************************
*
*     g4top.f:
*
*     This function returns the value of the inclusive structure function
*     g4.
*
************************************************************************
      function g4top(x)
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
      double precision g4top
*
      if(Polarized)then
*      
         if(TMC)then
            write(6,*) "TMCs not available for polarised DIS"
            call exit(-10)
         else
            if(x.lt.xmin(1)-tol.or.x.gt.xmax+tol)then
               write(6,*) "In g4top.f:"
               write(6,*) "Invalid value of x =",x
               call exit(-10)
            endif
            if (x.lt.xmin(1)) x = xmin(1)
            if (x.gt.xmax) x = 1d0
*     
*     Interpolation
*     
            g4top = 0d0
            n = inter_degree(0)
            do alpha=0,nin(0)
               g4top = g4top + w_int_gen(n,alpha,x) * F3(6,0,alpha)
            enddo
            if(dabs(g4top).le.1d-14) g4top = 0d0
         endif
      else
         write(6,*) "g4 structure function not available",
     1        " for unpolarised DIS"
         call exit(-10)
      endif
*
      return
      end
