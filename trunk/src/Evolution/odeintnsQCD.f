************************************************************************
*
*     odeintnsQCD.f:
*
*     Set of routines that perform the Adaptive Stepsize Control for 
*     Runge-Kutta integration of an Ordinary Differential Equation.
*
************************************************************************
      subroutine odeintnsQCD(i,mu21,mu22,ystart,y)
*
      implicit none
*
      include "../commons/PDFEvolution.h"
      include "../commons/grid.h"
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
      integer maxstp
      double precision x1,x2
      double precision a_QCD
      double precision h1,eps
      double precision h,hdid,hnext,x
      double precision dydx(0:nint_max,0:nint_max)
      double precision yscal(0:nint_max,0:nint_max)
      double precision tiny

      parameter(maxstp=1000)
      parameter(tiny=1d-10)
      parameter(h1=1d-3)
      parameter(eps=1d-3)
**
*     Output Variables
*
      double precision y(0:nint_max,0:nint_max)
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
         do beta=0,nin(igrid)
            y(alpha,beta) = ystart(alpha,beta)
         enddo
      enddo
*
      do nstp=1,maxstp
         call derivsnsQCD(i,x,y,dydx)
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
         call rkqsnsQCD(i,y,dydx,x,h,eps,yscal,hdid,hnext)
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
      subroutine rkqsnsQCD(i,y,dydx,x,htry,eps,yscal,hdid,hnext)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Variables
*
      integer i
      integer alpha,beta
      double precision safety,pgrow,pshrnk,errcon
      double precision eps,hdid,hnext,htry,x
      double precision dydx(0:nint_max,0:nint_max)
      double precision y(0:nint_max,0:nint_max)
      double precision yscal(0:nint_max,0:nint_max)
      double precision errmax,h,htemp,xnew
      double precision yerr(0:nint_max,0:nint_max)
      double precision ytemp(0:nint_max,0:nint_max)

      parameter(safety=0.9d0,pgrow=-0.2d0,pshrnk=-0.25d0,errcon=1.89d-4)
*
      h = htry
*
 101  call rkcknsQCD(i,y,dydx,x,h,ytemp,yerr)
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
      subroutine rkcknsQCD(i,y,dydx,x,h,yout,yerr)
*
      implicit none
*
      include "../commons/grid.h"
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
      double precision A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,
     1B51,B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,
     2DC6
      parameter(A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     1B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     2B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     3B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,
     4C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,
     5DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,
     6DC6=C6-.25)
*
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            ytemp(alpha,beta) = y(alpha,beta)
     1                        + B21 * h * dydx(alpha,beta)
         enddo
      enddo
      call derivsnsQCD(i,x+A2*h,ytemp,ak2)
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
           ytemp(alpha,beta) = y(alpha,beta) 
     1                       + h * ( B31 * dydx(alpha,beta)
     2                       +       B32 * ak2(alpha,beta) )
         enddo
      enddo
      call derivsnsQCD(i,x+A3*h,ytemp,ak3)
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            ytemp(alpha,beta) = y(alpha,beta) 
     1                        + h * ( B41 * dydx(alpha,beta)
     2                        +       B42 * ak2(alpha,beta) 
     3                        +       B43 * ak3(alpha,beta) )
         enddo
      enddo
      call derivsnsQCD(i,x+A4*h,ytemp,ak4)
      do alpha=0,nin(igrid)
         do beta=0,nin(igrid)
            ytemp(alpha,beta) = y(alpha,beta) 
     1                        + h * ( B51 * dydx(alpha,beta)
     2                        +       B52 * ak2(alpha,beta)
     3                        +       B53 * ak3(alpha,beta)
     4                        +       B54 * ak4(alpha,beta) )
         enddo
      enddo
      call derivsnsQCD(i,x+A5*h,ytemp,ak5)
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
      call derivsnsQCD(i,x+A6*h,ytemp,ak6)
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
      subroutine derivsnsQCD(i,t,M,dMdt)
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
      double precision M(0:nint_max,0:nint_max)
**
*     Internal Variables
*
      integer alpha,beta,delta
      double precision mu2
      double precision integralsQCD
      double precision coup,a_QCD
      double precision integ(0:nint_max)
**
*     Output Variables
*
      double precision dMdt(0:nint_max,0:nint_max)
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
*
************************************************************************
*
*     The following routine computes the derivative of the non-singlet
*     evolution operators in QCD.
*
************************************************************************
      subroutine DeriveNsQCD(i,coup,dMdt)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      integer i
      double precision coup
**
*     Internal Variables
*
      integer alpha,beta
      double precision integralsQCD
      double precision integ(0:nint_max)
**
*     Output Variables
*
      double precision dMdt(0:nint_max,0:nint_max)
*
      do alpha=0,nin(igrid)
         integ(alpha) = integralsQCD(0,alpha,coup,i)
      enddo
*
*     Initialization
*
      do alpha=0,nin(igrid)
         do beta=0,alpha-1
            dMdt(alpha,beta) = 0d0
         enddo
         do beta=alpha,nin(igrid)
            dMdt(alpha,beta) = integ(beta-alpha)
         enddo
      enddo
*
      return
      end
