************************************************************************
*
*     odeintsgQEDf.f:
*
*     Set of routines that perform the Adaptive Stepsize Control for 
*     Runge-Kutta integration of an Ordinary Differential Equation.
*
************************************************************************
      subroutine odeintsgQEDf(mu21,mu22,ystart,y)
*
      implicit none
*
      include "../commons/PDFEvolution.h"
      include "../commons/grid.h"
      include "../commons/odeint1.h"
**
*     Input Variables
*
      double precision mu21,mu22
      double precision ystart(3,0:nint_max)
**
*     Internal Variables
*
      integer i,nstp
      integer alpha
      double precision x1,x2
      double precision a_QED
      double precision h,hdid,hnext,x
      double precision dydx(3,0:nint_max)
      double precision yscal(3,0:nint_max)
**
*     Output Variables
*
      double precision y(3,0:nint_max)
*
      if(PDFEvol.eq."exactmu")then
         x1 = dlog(mu21)
         x2 = dlog(mu22)
      else
         x1 = a_QED(mu21)
         x2 = a_QED(mu22)
      endif
*
      x = x1
      h = sign(h1,x2-x1)
*
      do i=1,3
         do alpha=0,nin(igrid)
            y(i,alpha) = ystart(i,alpha)
         enddo
      enddo
*
      do nstp=1,maxstp
         call derivssgQEDf(x,y,dydx)
*
         do i=1,3
            do alpha=0,nin(igrid)
               yscal(i,alpha) = dabs(y(i,alpha)) 
     1                        + dabs(h*dydx(i,alpha)) 
     2                        + tiny
            enddo
         enddo
*
         if((x+h-x2)*(x+h-x1).gt.0d0) h = x2 - x
*
         call rkqssgQEDf(y,dydx,x,h,eps,yscal,hdid,hnext)
*
         if((x-x2)*(x2-x1).ge.0d0) return
*
         h = hnext
      enddo
*
      write(6,*) "In odeintsg.f:"
      write(6,*) "too many steps!"
      call exit(-10)
*
      return
      end
*
************************************************************************
      subroutine rkqssgQEDf(y,dydx,x,htry,eps,yscal,hdid,hnext)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/odeint2.h"
**
*     Variables
*
      integer i
      integer alpha
      double precision eps,hdid,hnext,htry,x
      double precision dydx(3,0:nint_max)
      double precision y(3,0:nint_max)
      double precision yscal(3,0:nint_max)
      double precision errmax,h,htemp,xnew
      double precision yerr(3,0:nint_max)
      double precision ytemp(3,0:nint_max)
*
      h = htry
*
 101  call rkcksgQEDf(y,dydx,x,h,ytemp,yerr)
*
      errmax = 0d0
      do i=1,3
         do alpha=0,nin(igrid)
            errmax = max(errmax,dabs(yerr(i,alpha)
     1                  /yscal(i,alpha)))
         enddo
      enddo
*
      errmax = errmax / eps
*
      if(errmax.gt.1d0)then
         htemp = safety * h * (errmax**pshrnk)
         h     = sign(max(dabs(htemp),0.1d0*dabs(h)),h)
         xnew  = x + h
         if(xnew.eq.x)then
            write(6,*) "In odeintsg.f:"
            write(6,*) "stepsize underflow in rkqssg"
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
         do i=1,3
            do alpha=0,nin(igrid)
               y(i,alpha) = ytemp(i,alpha)
            enddo
         enddo
         return
      endif
*
      end
*
************************************************************************
      subroutine rkcksgQEDf(y,dydx,x,h,yout,yerr)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/odeint2.h"
**
*     Variables
*
      integer i
      integer alpha
      double precision h,x
      double precision dydx(3,0:nint_max)
      double precision y(3,0:nint_max)
      double precision yerr(3,0:nint_max)
      double precision yout(3,0:nint_max)
      double precision ytemp(3,0:nint_max)
      double precision ak2(3,0:nint_max)
      double precision ak3(3,0:nint_max)
      double precision ak4(3,0:nint_max)
      double precision ak5(3,0:nint_max)
      double precision ak6(3,0:nint_max)
*
      do i=1,3
         do alpha=0,nin(igrid)
            ytemp(i,alpha) = y(i,alpha)
     1                     + B21 * h * dydx(i,alpha)
         enddo
      enddo
      call derivssgQEDf(x+A2*h,ytemp,ak2)
      do i=1,3
         do alpha=0,nin(igrid)
            ytemp(i,alpha) = y(i,alpha) 
     1                     + h * ( B31 * dydx(i,alpha)
     2                     +       B32 * ak2(i,alpha) )
         enddo
      enddo
      call derivssgQEDf(x+A3*h,ytemp,ak3)
      do i=1,3
         do alpha=0,nin(igrid)
            ytemp(i,alpha) = y(i,alpha) 
     1                     + h * ( B41 * dydx(i,alpha)
     2                     +       B42 * ak2(i,alpha) 
     3                     +       B43 * ak3(i,alpha) )
         enddo
      enddo
      call derivssgQEDf(x+A4*h,ytemp,ak4)
      do i=1,3
         do alpha=0,nin(igrid)
            ytemp(i,alpha) = y(i,alpha) 
     1                     + h * ( B51 * dydx(i,alpha)
     2                     +       B52 * ak2(i,alpha)
     3                     +       B53 * ak3(i,alpha)
     4                     +       B54 * ak4(i,alpha) )
         enddo
      enddo
      call derivssgQEDf(x+A5*h,ytemp,ak5)
      do i=1,3
         do alpha=0,nin(igrid)
            ytemp(i,alpha) = y(i,alpha)
     1                     + h * ( B61 * dydx(i,alpha)
     2                     +       B62 * ak2(i,alpha)
     3                     +       B63 * ak3(i,alpha)
     4                     +       B64 * ak4(i,alpha)
     5                     +       B65 * ak5(i,alpha) )
         enddo
      enddo
      call derivssgQEDf(x+A6*h,ytemp,ak6)
      do i=1,3
         do alpha=0,nin(igrid)
            yout(i,alpha) = y(i,alpha)
     1                    + h * ( C1 * dydx(i,alpha)
     2                    +       C3 * ak3(i,alpha)
     3                    +       C4 * ak4(i,alpha)
     4                    +       C6 * ak6(i,alpha) )
         enddo
      enddo
      do i=1,3
         do alpha=0,nin(igrid)
            yerr(i,alpha) = h *( DC1 * dydx(i,alpha)
     1                    +      DC3 * ak3(i,alpha)
     2                    +      DC4 * ak4(i,alpha)
     3                    +      DC5 * ak5(i,alpha)
     4                    +      DC6 * ak6(i,alpha) )
         enddo
      enddo
*
      return
      end
*
************************************************************************
      subroutine derivssgQEDf(t,f,dfdt)
*
      implicit none
*
      include "../commons/PDFEvolution.h"
      include "../commons/grid.h"
**
*     Input Variables
*
      double precision t
      double precision f(3,0:nint_max)
**
*     Internal Variables
*
      integer i,k,mapp(3,3)
      integer alpha,delta
      double precision mu2
      double precision integralsQED
      double precision coup,a_QED
      double precision integ(0:nint_max,3,3)
**
*     Output Variables
*
      double precision dfdt(3,0:nint_max)
*
      if(PDFEvol.eq."exactmu")then
         mu2  = dexp(t)
         coup = a_QED(mu2)
      else
         coup = t
      endif
*
*     Map used for the muliplication
*
      mapp(1,1) = 3
      mapp(1,2) = 4
      mapp(1,3) = 5
      mapp(2,1) = 6
      mapp(2,2) = 7
      mapp(2,3) = 8
      mapp(3,1) = 9
      mapp(3,2) = 10
      mapp(3,3) = 11
*
      do alpha=0,nin(igrid)
         do i=1,3
            do k=1,3
               integ(alpha,i,k) =
     1              integralsQED(0,alpha,coup,0d0,mapp(i,k))
            enddo
         enddo
      enddo
*
*     Initialization
*
      do i=1,3
         do alpha=0,nin(igrid)
            dfdt(i,alpha) = 0d0
            do k=1,3
               do delta=0,nin(igrid)-alpha
                  dfdt(i,alpha) = dfdt(i,alpha)
     1            + integ(delta,i,k) * f(k,alpha+delta)
               enddo
            enddo
         enddo
      enddo
*
      return
      end
