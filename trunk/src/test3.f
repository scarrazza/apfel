************************************************************************
*
*     program test3.f:
*
*     test of the Runge-Kutta routine.
*
************************************************************************
      program test3
*
      implicit none
*
*     Variables
*
      integer nvar
      integer i,j
      double precision h1,hmin,eps
      double precision x1,x2
      double precision ystart(10)
      double precision y(10)
      double precision step
      external polynomial
      external exponential
      external coupledexponential
*
*     Equation:
*
*      dy
*     ---- = x
*      dx
*
      nvar = 1
      h1   = 0.001d0
      hmin = 0d0
      eps  = 1d-5
      x1   = 0d0
      x2   = x1
*
*     Initial conditions
*
      do i=1,nvar
         ystart(i) = 0d0
      enddo
      step = 1d-3
      do j=1,10001
         call odeint(nvar,h1,hmin,eps,x1,x2,ystart,y,polynomial)
*
*        expected solution y = x^2 / 2
*
         write(6,*) x1,x2,(dabs(y(i)-x2**2d0/2d0),i=1,nvar)
         x2 = x2 + step
      enddo
*
      write(6,*) "****************************"
*
*     Equation:
*
*      dy
*     ---- = y
*      dx
*
      nvar = 1
      h1   = 0.0001d0
      hmin = 0d0
      eps  = 1d-5
      x1   = 0d0
      x2   = x1
*
*     Initial conditions
*
      do i=1,nvar
         ystart(i) = 1d0
      enddo
      step = 1d-3
      do j=1,1001
         call odeint(nvar,h1,hmin,eps,x1,x2,ystart,y,exponential)
*
*     Expected solution y = exp(x)
*
         write(6,*) x1,x2,(dabs(y(i)-dexp(x2)),i=1,nvar)
         x2 = x2 + step
      enddo
*
      write(6,*) "****************************"
*
*     Equation:
*      _
*     |  dy(1)
*     |  ---- = y(2)
*     |   dx
*     |
*     |  dy(2)
*     |  ---- = y(1)
*     |_  dx
*
      nvar = 2
      h1   = 0.0001d0
      hmin = 0d0
      eps  = 1d-5
      x1   = 0d0
      x2   = x1
*
*     Initial conditions
*
      do i=1,nvar
         ystart(i) = 1d0
      enddo
      step = 1d-3
      do j=1,1001
         call odeint(nvar,h1,hmin,eps,x1,x2,ystart,y,coupledexponential)
*
*     Expected solution y(1) = y(2) = exp(x)
*
         write(6,*) x1,x2,(dabs(y(i)-dexp(x2)),i=1,nvar)
         x2 = x2 + step
      enddo
*
      return
      end
