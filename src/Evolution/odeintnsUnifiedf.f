************************************************************************
*
*     odeintnsUnifiedf.f:
*
*     Set of routines that perform the Adaptive Stepsize Control for 
*     Runge-Kutta integration of an Ordinary Differential Equation.
*
************************************************************************
      subroutine odeintnsUnifiedf(i,mu21,mu22,ystart,y)
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
         call derivsnsUnifiedf(i,x,y,dydx)
*
         do alpha=0,nin(igrid)
            yscal(alpha) = dabs(y(alpha)) + dabs(h*dydx(alpha)) + tiny
         enddo
*
         if((x+h-x2)*(x+h-x1).gt.0d0) h = x2 - x
*
         call rkqsnsUnifiedf(i,y,dydx,x,h,eps,yscal,hdid,hnext)
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
      subroutine rkqsnsUnifiedf(i,y,dydx,x,htry,eps,yscal,hdid,hnext)
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
 101  call rkcknsUnifiedf(i,y,dydx,x,h,ytemp,yerr)
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
      subroutine rkcknsUnifiedf(i,y,dydx,x,h,yout,yerr)
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
      call derivsnsUnifiedf(i,x+A2*h,ytemp,ak2)
      do alpha=0,nin(igrid)
         ytemp(alpha) = y(alpha) 
     1                + h * ( B31 * dydx(alpha)
     2                +       B32 * ak2(alpha) )
      enddo
      call derivsnsUnifiedf(i,x+A3*h,ytemp,ak3)
      do alpha=0,nin(igrid)
         ytemp(alpha) = y(alpha) 
     1                + h * ( B41 * dydx(alpha)
     2                +       B42 * ak2(alpha) 
     3                +       B43 * ak3(alpha) )
      enddo
      call derivsnsUnifiedf(i,x+A4*h,ytemp,ak4)
      do alpha=0,nin(igrid)
         ytemp(alpha) = y(alpha) 
     1                + h * ( B51 * dydx(alpha)
     2                +       B52 * ak2(alpha)
     3                +       B53 * ak3(alpha)
     4                +       B54 * ak4(alpha) )
      enddo
      call derivsnsUnifiedf(i,x+A5*h,ytemp,ak5)
      do alpha=0,nin(igrid)
         ytemp(alpha) = y(alpha)
     1                + h * ( B61 * dydx(alpha)
     2                +       B62 * ak2(alpha)
     3                +       B63 * ak3(alpha)
     4                +       B64 * ak4(alpha)
     5                +       B65 * ak5(alpha) )
      enddo
      call derivsnsUnifiedf(i,x+A6*h,ytemp,ak6)
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
      subroutine derivsnsUnifiedf(i,t,f,dfdt)
*
      implicit none
*
      include "../commons/PDFEvolution.h"
      include "../commons/grid.h"
      include "../commons/wrap.h"
      include "../commons/ipt.h"
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
      double precision integralsQED
      double precision coupQCD,a_QCD,muR2,bts,fbeta
      double precision coupQED,a_QED
      double precision integ1(0:nint_max)
      double precision integ2(0:nint_max,0:nint_max)
**
*     Output Variables
*
      double precision dfdt(0:nint_max)
*
      if(PDFEvol.eq."exactmu")then
         mu2     = dexp(t)
         coupQCD = a_QCD(mu2)
         coupQED = a_QED(mu2)
         bts     = 1d0
      else
         mu2     = muR2(t)
         coupQCD = t
         coupQED = a_QED(mu2)
         bts     = 1d0 / fbeta(t,wnf,ipt)
      endif
*
      if(IsExt(igrid))then
         if(i.eq.1)then
            do alpha=0,nin(igrid)
               do delta=alpha,nin(igrid)
                  integ2(alpha,delta) =
     1                 integralsQCD(alpha,delta,coupQCD,1) + bts
     2                 * integralsQED(alpha,delta,coupQED,coupQCD,1)
               enddo
            enddo
         elseif(i.eq.2)then
            do alpha=0,nin(igrid)
               do delta=alpha,nin(igrid)
                  integ2(alpha,delta) =
     1                 integralsQCD(alpha,delta,coupQCD,1) + bts
     2                 * integralsQED(alpha,delta,coupQED,coupQCD,2)
               enddo
            enddo
         elseif(i.eq.3)then
            do alpha=0,nin(igrid)
               do delta=alpha,nin(igrid)
                  integ2(alpha,delta) =
     1                 integralsQCD(alpha,delta,coupQCD,2) + bts
     2                 * integralsQED(alpha,delta,coupQED,coupQCD,3)
               enddo
            enddo
         elseif(i.eq.4)then
            do alpha=0,nin(igrid)
               do delta=alpha,nin(igrid)
                  integ2(alpha,delta) =
     1                 integralsQCD(alpha,delta,coupQCD,2) + bts
     2                 * integralsQED(alpha,delta,coupQED,coupQCD,4)
               enddo
            enddo
         elseif(i.eq.5)then
            do alpha=0,nin(igrid)
               do delta=alpha,nin(igrid)
               integ2(alpha,delta) = bts
     1                 * integralsQED(alpha,delta,coupQED,coupQCD,23)
               enddo
            enddo
         endif
*
*     Convolution
*
         do alpha=0,nin(igrid)
            dfdt(alpha) = 0d0
            do delta=alpha,nin(igrid)
               dfdt(alpha) = dfdt(alpha)
     1              + integ2(alpha,delta) * f(delta)
            enddo
         enddo
      else
         if(i.eq.1)then
            do alpha=0,nin(igrid)
               integ1(alpha) = integralsQCD(0,alpha,coupQCD,1)
     1              + bts * integralsQED(0,alpha,coupQED,coupQCD,1)
            enddo
         elseif(i.eq.2)then
            do alpha=0,nin(igrid)
               integ1(alpha) = integralsQCD(0,alpha,coupQCD,1)
     1              + bts * integralsQED(0,alpha,coupQED,coupQCD,2)
            enddo
         elseif(i.eq.3)then
            do alpha=0,nin(igrid)
               integ1(alpha) = integralsQCD(0,alpha,coupQCD,2)
     1              + bts * integralsQED(0,alpha,coupQED,coupQCD,3)
            enddo
         elseif(i.eq.4)then
            do alpha=0,nin(igrid)
               integ1(alpha) = integralsQCD(0,alpha,coupQCD,2)
     1              + bts * integralsQED(0,alpha,coupQED,coupQCD,4)
            enddo
         elseif(i.eq.5)then
            do alpha=0,nin(igrid)
               integ1(alpha) = bts
     1              * integralsQED(0,alpha,coupQED,coupQCD,23)
            enddo
         endif
*
*     Convolution
*
         do alpha=0,nin(igrid)
            dfdt(alpha) = 0d0
            do delta=0,nin(igrid)-alpha
               dfdt(alpha) = dfdt(alpha)
     1              + integ1(delta) * f(alpha+delta)
            enddo
         enddo
      endif
*
      return
      end
