************************************************************************
*
*     odeintnsUnified.f:
*
*     Set of routines that perform the Adaptive Stepsize Control for 
*     Runge-Kutta integration of an Ordinary Differential Equation.
*
************************************************************************
      subroutine odeintnsUnified(i,mu21,mu22,ystart,y)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/odeint1.h"
**
*     Input Variables
*
      integer i
      double precision mu21,mu22
      double precision ystart(0:nint_max,0:nint_max)
**
*     Internal Variables
*
      integer nstp
      integer alpha,beta
      double precision x1,x2
      double precision h,hdid,hnext,x
      double precision dydx(0:nint_max,0:nint_max)
      double precision yscal(0:nint_max,0:nint_max)
**
*     Output Variables
*
      double precision y(0:nint_max,0:nint_max)
*
      x1 = dlog(mu21)
      x2 = dlog(mu22)
*
      x = x1
      h = sign(h1,x2-x1)
*
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            y(alpha,beta) = ystart(alpha,beta)
         enddo
      enddo
*
      do nstp=1,maxstp
         call derivsnsUnified(i,x,y,dydx)
*
         do alpha=0,nin(igrid)
            do beta=0,nin(igrid)
               yscal(alpha,beta) = dabs(y(alpha,beta)) 
     1                           + dabs(h*dydx(alpha,beta)) 
     2                           + tiny
            enddo
         enddo
*
         if((x+h-x2)*(x+h-x1).gt.0d0) h = x2 - x
*
         call rkqsnsUnified(i,y,dydx,x,h,eps,yscal,hdid,hnext)
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
      subroutine rkqsnsUnified(i,y,dydx,x,htry,eps,yscal,hdid,hnext)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/odeint2.h"
**
*     Variables
*
      integer i
      integer alpha,beta
      double precision eps,hdid,hnext,htry,x
      double precision dydx(0:nint_max,0:nint_max)
      double precision y(0:nint_max,0:nint_max)
      double precision yscal(0:nint_max,0:nint_max)
      double precision errmax,h,htemp,xnew
      double precision yerr(0:nint_max,0:nint_max)
      double precision ytemp(0:nint_max,0:nint_max)
*
      h = htry
*
 101  call rkcknsUnified(i,y,dydx,x,h,ytemp,yerr)
*
      errmax = 0d0
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            errmax = max(errmax,dabs(yerr(alpha,beta)
     1             /yscal(alpha,beta)))
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
            do beta=0,nin(igrid)
               y(alpha,beta) = ytemp(alpha,beta)
            enddo
         enddo
         return
      endif
*
      end
*
************************************************************************
      subroutine rkcknsUnified(i,y,dydx,x,h,yout,yerr)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/odeint2.h"
**
*     Variables
*
      integer i
      integer alpha,beta
      double precision h,x
      double precision dydx(0:nint_max,0:nint_max)
      double precision y(0:nint_max,0:nint_max)
      double precision yerr(0:nint_max,0:nint_max)
      double precision yout(0:nint_max,0:nint_max)
      double precision ytemp(0:nint_max,0:nint_max)
      double precision ak2(0:nint_max,0:nint_max)
      double precision ak3(0:nint_max,0:nint_max)
      double precision ak4(0:nint_max,0:nint_max)
      double precision ak5(0:nint_max,0:nint_max)
      double precision ak6(0:nint_max,0:nint_max)
*
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            ytemp(alpha,beta) = y(alpha,beta)
     1                        + B21 * h * dydx(alpha,beta)
         enddo
      enddo
      call derivsnsUnified(i,x+A2*h,ytemp,ak2)
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
           ytemp(alpha,beta) = y(alpha,beta) 
     1                       + h * ( B31 * dydx(alpha,beta)
     2                       +       B32 * ak2(alpha,beta) )
         enddo
      enddo
      call derivsnsUnified(i,x+A3*h,ytemp,ak3)
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            ytemp(alpha,beta) = y(alpha,beta) 
     1                        + h * ( B41 * dydx(alpha,beta)
     2                        +       B42 * ak2(alpha,beta) 
     3                        +       B43 * ak3(alpha,beta) )
         enddo
      enddo
      call derivsnsUnified(i,x+A4*h,ytemp,ak4)
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            ytemp(alpha,beta) = y(alpha,beta) 
     1                        + h * ( B51 * dydx(alpha,beta)
     2                        +       B52 * ak2(alpha,beta)
     3                        +       B53 * ak3(alpha,beta)
     4                        +       B54 * ak4(alpha,beta) )
         enddo
      enddo
      call derivsnsUnified(i,x+A5*h,ytemp,ak5)
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            ytemp(alpha,beta) = y(alpha,beta)
     1                        + h * ( B61 * dydx(alpha,beta)
     2                        +       B62 * ak2(alpha,beta)
     3                        +       B63 * ak3(alpha,beta)
     4                        +       B64 * ak4(alpha,beta)
     5                        +       B65 * ak5(alpha,beta) )
         enddo
      enddo
      call derivsnsUnified(i,x+A6*h,ytemp,ak6)
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            yout(alpha,beta) = y(alpha,beta)
     1                       + h * ( C1 * dydx(alpha,beta)
     2                       +       C3 * ak3(alpha,beta)
     3                       +       C4 * ak4(alpha,beta)
     4                       +       C6 * ak6(alpha,beta) )
         enddo
      enddo
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            yerr(alpha,beta) = h * ( DC1 * dydx(alpha,beta)
     1                       +       DC3 * ak3(alpha,beta)
     2                       +       DC4 * ak4(alpha,beta)
     3                       +       DC5 * ak5(alpha,beta)
     4                       +       DC6 * ak6(alpha,beta) )
         enddo
      enddo
*
      return
      end
*
************************************************************************
      subroutine derivsnsUnified(i,t,M,dMdt)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      integer i
      double precision t
      double precision M(0:nint_max,0:nint_max)
**
*     Internal Variables
*
      integer alpha,beta,delta
      double precision mu2
      double precision integralsQCD
      double precision integralsQED
      double precision coupQCD,a_QCD
      double precision coupQED,a_QED
      double precision integ(0:nint_max)
**
*     Output Variables
*
      double precision dMdt(0:nint_max,0:nint_max)
*
      mu2     = dexp(t)
      coupQCD = a_QCD(mu2)
      coupQED = a_QED(mu2)
*
      if(i.eq.1)then
         do alpha=0,nin(igrid)
            integ(alpha) = integralsQCD(0,alpha,coupQCD,1)
     1                   + integralsQED(0,alpha,coupQED,1)
         enddo
      elseif(i.eq.2)then
         do alpha=0,nin(igrid)
            integ(alpha) = integralsQCD(0,alpha,coupQCD,1)
     1                   + integralsQED(0,alpha,coupQED,2)
         enddo
      elseif(i.eq.3)then
         do alpha=0,nin(igrid)
            integ(alpha) = integralsQCD(0,alpha,coupQCD,2)
     1                   + integralsQED(0,alpha,coupQED,1)
         enddo
      elseif(i.eq.4)then
         do alpha=0,nin(igrid)
            integ(alpha) = integralsQCD(0,alpha,coupQCD,2)
     1                   + integralsQED(0,alpha,coupQED,2)
         enddo
      endif
*
*     Initialization
*
      do alpha=0,nin(igrid)
         do beta=alpha,nin(igrid)
            dMdt(alpha,beta) = 0d0
            do delta=0,nin(igrid)-alpha
               dMdt(alpha,beta) = dMdt(alpha,beta)
     1                          + integ(delta) * M(alpha+delta,beta)
            enddo
         enddo
      enddo
*
      return
      end
