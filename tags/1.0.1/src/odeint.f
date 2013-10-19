************************************************************************
*
*     odeint.f:
*
*     Set of routines that perform the Adaptive Stepsize Control for 
*     Runge-Kutta integration of an Ordinary Differential Equation.
*
*     Integrate the starting values ystart(1:nvar) from x1 to x2 with 
*     accuracy eps.
*
*     h1 should be set as a guessed first stepsize, hmin as the minimum 
*     allowed stepsize (can be zero).
*
*     derivs is the user-supplied subroutine for calculating the 
*     right-hand side derivative, while rkqs is the name of the stepper 
*     routine to be used.
*
************************************************************************
      subroutine odeint(nvar,h1,hmin,eps,x1,x2,ystart,y,derivs)
*
      implicit none
*
      include "../commons/rk.h"
**
*     Input Variables
*
      integer nvar
      double precision h1,hmin,eps
      double precision x1,x2
      double precision ystart(nvar)
      external derivs
**
*     Internal Variables
*
      integer i,nstp
      double precision h,hdid,hnext,x
      double precision dydx(nmax)
      double precision yscal(nmax)
**
*     Output Variables
*
      double precision y(nvar)
*
      x = x1
      h = sign(h1,x2-x1)
*
      do i=1,nvar
         y(i) = ystart(i)
      enddo
*
      do nstp=1,maxstp
         call derivs(nvar,x,y,dydx)
*
         do i=1,nvar
            yscal(i) = abs(y(i)) + abs( h * dydx(i) ) + tiny
         enddo
*
         if((x+h-x2)*(x+h-x1).gt.0d0) h = x2 - x
*
         call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
*
         if((x-x2)*(x2-x1).ge.0d0) return
*
         if(abs(hnext).lt.hmin)then
            write(6,*) "In odeint.f:"
            write(6,*) "stepsize smaller than minimum in odeint"
            call exit(-10)
         endif
         h = hnext
      enddo
*
      write(6,*) "In odeint.f:"
      write(6,*) "too many steps!"
      call exit(-10)
*
      return
      end
*
************************************************************************
      subroutine rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
*
      implicit none
*
      include "../commons/rk.h"
**
*     Variables
*
      integer n
      integer i
      double precision eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      double precision errmax,h,htemp,xnew,yerr(nmax),ytemp(nmax)
      external derivs
*
      h = htry
*
 101  call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
*
      errmax = 0d0
      do i=1,n
         errmax = max(errmax,abs(yerr(i)/yscal(i)))
      enddo
*
      errmax = errmax / eps
*
      if(errmax.gt.1d0)then
         htemp = safety * h * (errmax**pshrnk)
         h     = sign(max(abs(htemp),0.1d0*abs(h)),h)
         xnew  = x + h
         if(xnew.eq.x)then
            write(6,*) "In odeint.f:"
            write(6,*) "stepsize underflow in rkqs"
            call exit(-10)
         endif
         goto 101
      else
         if(errmax.gt.errcon)then
            hnext = safety * h * (errmax**pgrow)
         else
            hnext = 5d0 * h
         endif
         hdid = h
         x    = x + h
         do i=1,n
            y(i) = ytemp(i)
         enddo
         return
      endif
*
      end
*
************************************************************************
      subroutine rkck(y,dydx,n,x,h,yout,yerr,derivs)
*
      implicit none
*
      include "../commons/rk.h"
**
*     Variables
*
      integer n
      integer i
      double precision h,x,dydx(n),y(n),yerr(n),yout(n)
      double precision ak2(nmax),ak3(nmax),ak4(nmax),ak5(nmax),
     1ak6(nmax),ytemp(nmax),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,
     2B51,B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,
     3DC6
      parameter(A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     1B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     2B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     3B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,
     4C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,
     5DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,
     6DC6=C6-.25)
      external derivs
*
      do 11 i=1,n
         ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
      call derivs(n,x+A2*h,ytemp,ak2)
      do 12 i=1,n
         ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
      call derivs(n,x+A3*h,ytemp,ak3)
      do 13 i=1,n
         ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(n,x+A4*h,ytemp,ak4)
      do 14 i=1,n
         ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(n,x+A5*h,ytemp,ak5)
      do 15 i=1,n
         ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     1            B65*ak5(i))
15    continue
      call derivs(n,x+A6*h,ytemp,ak6)
      do 16 i=1,n
         yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
         yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     1           ak6(i))
17    continue
*
      return
      end
