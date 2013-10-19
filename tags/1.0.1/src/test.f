************************************************************************
*
*      test program to check the interpolant functions.
*
************************************************************************
      program test
*
      implicit none
*
      include "../commons/grid.h"
*
      integer n,alpha
      double precision w_int,x,x_int
*
*     Read input parameters
*
      call initParameters
*
*     Initialize grid
*
      call initGrid
*
*     Test of the interpolant routine
*
      write(6,*) "Enter the value of x"
      read(5,*) x
*
      n = inter_degree
*
      x_int = 0d0
      do alpha=0,nint+inter_degree
         x_int = x_int + w_int(n,alpha,x) * xg(alpha)
c         write(6,*) alpha,xg(alpha),w_int(n,alpha,x)
      enddo
      write(6,*) "  "
*
      write(6,*) "       Input       ",
     1           "    Intepolated    ",
     2           "     Accuracy      "
      write(6,'(E18.10,1X,E18.10,1X,E18.10)') x,x_int,dabs(x-x_int)/x
      write(6,*) "  "
*
      end
