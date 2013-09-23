************************************************************************
*
*     coupledexponential.f:
*
*     It returns the r.h.s. of the differential equation:
*      _
*     |  dy(1)
*     |  ---- = y(2)
*     |   dx
*     |
*     |  dy(2)
*     |  ---- = y(1)
*     |_  dx
*
*     to be used in the program test3.f to the Runge-Kutta method.
*
************************************************************************
      subroutine coupledexponential(n,x,y,dydx)
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
         dydx(i) = y(n+1-i)
      enddo
*
      return
      end
