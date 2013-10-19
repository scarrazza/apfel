************************************************************************
*
*     exponential.f:
*
*     It returns the r.h.s. of the differential equation:
*
*      dy
*     ---- = y
*      dx
*
*     to be used in the program test3.f to the Runge-Kutta method.
*
************************************************************************
      subroutine exponential(n,x,y,dydx)
*
      implicit none
**
*     Input Variables
*
      integer n
      double precision x
      double precision y(n)
**
*     Internal Variables
*
      integer i
**
*     Output Variables
*
      double precision dydx(n)
*
      do i=1,n
         dydx(i) = y(i)
      enddo
*
      return
      end
