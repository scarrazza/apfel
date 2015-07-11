************************************************************************
*
*     odeintnsQCDf.f:
*
*     Set of routines that perform the Adaptive Stepsize Control for 
*     Runge-Kutta integration of an Ordinary Differential Equation.
*
************************************************************************
      subroutine odeintnsQCDf(i,mu21,mu22,ystart,y)
*
      implicit none
*
      include "../commons/PDFEvolution.h"
      include "../commons/grid.h"
      include "../commons/odeint1.h"
**
*     Input Variables
*
      integer i
      double precision mu21,mu22
      double precision ystart(0:nint_max)
**
*     Internal Variables
*
      integer nstp
      integer alpha
      double precision x1,x2
      double precision a_QCD
      double precision h,hdid,hnext,x
      double precision dydx(0:nint_max)
      double precision yscal(0:nint_max)
**
*     Output Variables
*
      double precision y(0:nint_max)
*
      if(PDFEvol.eq."exactmu")then
         x1 = dlog(mu21)
         x2 = dlog(mu22)
      else
         x1 = a_QCD(mu21)
         x2 = a_QCD(mu22)
      endif
*
      x = x1
      h = sign(h1,x2-x1)
*
      do alpha=0,nin(igrid)
         y(alpha) = ystart(alpha)
      enddo
*
      do nstp=1,maxstp
         call derivsnsQCDf(i,x,y,dydx)
*
         do alpha=0,nin(igrid)
            yscal(alpha) = dabs(y(alpha)) + dabs(h*dydx(alpha)) + tiny
         enddo
*
         if((x+h-x2)*(x+h-x1).gt.0d0) h = x2 - x
*
         call rkqsnsQCDf(i,y,dydx,x,h,eps,yscal,hdid,hnext)
*
         if((x-x2)*(x2-x1).ge.0d0) return
*
         h = hnext
      enddo
*
      write(6,*) "In odeintns.f:"
      write(6,*) "too many steps!"
      call exit(-10)
*
      return
      end
*
************************************************************************
      subroutine rkqsnsQCDf(i,y,dydx,x,htry,eps,yscal,hdid,hnext)
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
      double precision dydx(0:nint_max)
      double precision y(0:nint_max)
      double precision yscal(0:nint_max)
      double precision errmax,h,htemp,xnew
      double precision yerr(0:nint_max)
      double precision ytemp(0:nint_max)
*
      h = htry
*
 101  call rkcknsQCDf(i,y,dydx,x,h,ytemp,yerr)
*
      errmax = 0d0
      do alpha=0,nin(igrid)
         errmax = max(errmax,dabs(yerr(alpha)/yscal(alpha)))
      enddo
*
      errmax = errmax / eps
*
      if(errmax.gt.1d0)then
         htemp = safety * h * (errmax**pshrnk)
         h     = sign(max(dabs(htemp),0.1d0*dabs(h)),h)
         xnew  = x + h
         if(xnew.eq.x)then
            write(6,*) "In odeintns.f:"
            write(6,*) "stepsize underflow in rkqsns"
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
         do alpha=0,nin(igrid)
            y(alpha) = ytemp(alpha)
         enddo
         return
      endif
*
      end
*
************************************************************************
      subroutine rkcknsQCDf(i,y,dydx,x,h,yout,yerr)
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
      double precision dydx(0:nint_max)
      double precision y(0:nint_max)
      double precision yerr(0:nint_max)
      double precision yout(0:nint_max)
      double precision ytemp(0:nint_max)
      double precision ak2(0:nint_max)
      double precision ak3(0:nint_max)
      double precision ak4(0:nint_max)
      double precision ak5(0:nint_max)
      double precision ak6(0:nint_max)
*
      do alpha=0,nin(igrid)
         ytemp(alpha) = y(alpha)
     1                + B21 * h * dydx(alpha)
      enddo
      call derivsnsQCDf(i,x+A2*h,ytemp,ak2)
      do alpha=0,nin(igrid)
         ytemp(alpha) = y(alpha) 
     1                + h * ( B31 * dydx(alpha)
     2                +       B32 * ak2(alpha) )
      enddo
      call derivsnsQCDf(i,x+A3*h,ytemp,ak3)
      do alpha=0,nin(igrid)
         ytemp(alpha) = y(alpha) 
     1                + h * ( B41 * dydx(alpha)
     2                +       B42 * ak2(alpha) 
     3                +       B43 * ak3(alpha) )
      enddo
      call derivsnsQCDf(i,x+A4*h,ytemp,ak4)
      do alpha=0,nin(igrid)
         ytemp(alpha) = y(alpha) 
     1                + h * ( B51 * dydx(alpha)
     2                +       B52 * ak2(alpha)
     3                +       B53 * ak3(alpha)
     4                +       B54 * ak4(alpha) )
      enddo
      call derivsnsQCDf(i,x+A5*h,ytemp,ak5)
      do alpha=0,nin(igrid)
         ytemp(alpha) = y(alpha)
     1                + h * ( B61 * dydx(alpha)
     2                +       B62 * ak2(alpha)
     3                +       B63 * ak3(alpha)
     4                +       B64 * ak4(alpha)
     5                +       B65 * ak5(alpha) )
      enddo
      call derivsnsQCDf(i,x+A6*h,ytemp,ak6)
      do alpha=0,nin(igrid)
         yout(alpha) = y(alpha)
     1               + h * ( C1 * dydx(alpha)
     2               +       C3 * ak3(alpha)
     3               +       C4 * ak4(alpha)
     4               +       C6 * ak6(alpha) )
      enddo
      do alpha=0,nin(igrid)
         yerr(alpha) = h * ( DC1 * dydx(alpha)
     1               +       DC3 * ak3(alpha)
     2               +       DC4 * ak4(alpha)
     3               +       DC5 * ak5(alpha)
     4               +       DC6 * ak6(alpha) )
      enddo
*
      return
      end
*
************************************************************************
      subroutine derivsnsQCDf(i,t,f,dfdt)
*
      implicit none
*
      include "../commons/PDFEvolution.h"
      include "../commons/grid.h"
**
*     Input Variables
*
      integer i
      double precision t
      double precision f(0:nint_max)
**
*     Internal Variables
*
      integer alpha,delta
      double precision mu2
      double precision integralsQCD
      double precision coup,a_QCD
      double precision integ(0:nint_max)
**
*     Output Variables
*
      double precision dfdt(0:nint_max)
*
      if(PDFEvol.eq."exactmu")then
         mu2  = dexp(t)
         coup = a_QCD(mu2)
      else
         coup = t
      endif
*
      do alpha=0,nin(igrid)
         integ(alpha) = integralsQCD(0,alpha,coup,i)
      enddo
*
*     Initialization
*
      do alpha=0,nin(igrid)
         dfdt(alpha) = 0d0
         do delta=0,nin(igrid)-alpha
            dfdt(alpha) = dfdt(alpha)
     1                  + integ(delta) * f(alpha+delta) 
         enddo
      enddo
*
      return
      end
