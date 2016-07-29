************************************************************************
*
*     odeintsgUnifiedf.f:
*
*     Set of routines that perform the Adaptive Stepsize Control for 
*     Runge-Kutta integration of an Ordinary Differential Equation.
*
************************************************************************
*
*     First Singolet
*
************************************************************************
      subroutine odeintsgUnifiedfS1(mu21,mu22,ystart,y)
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
      double precision ystart(5,0:nint_max)
**
*     Internal Variables
*
      integer i,nstp
      integer alpha
      double precision x1,x2
      double precision a_QCD
      double precision h,hdid,hnext,x
      double precision dydx(5,0:nint_max)
      double precision yscal(5,0:nint_max)
**
*     Output Variables
*
      double precision y(5,0:nint_max)
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
      do i=1,5
         do alpha=0,nin(igrid)
            y(i,alpha) = ystart(i,alpha)
         enddo
      enddo
*
      do nstp=1,maxstp
         call derivssgUnifiedfS1(x,y,dydx)
*
         do i=1,5
            do alpha=0,nin(igrid)
               yscal(i,alpha) = dabs(y(i,alpha)) 
     1                        + dabs(h*dydx(i,alpha)) 
     2                        + tiny
            enddo
         enddo
*
         if((x+h-x2)*(x+h-x1).gt.0d0) h = x2 - x
*
         call rkqssgUnifiedfS1(y,dydx,x,h,eps,yscal,hdid,hnext)
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
      subroutine rkqssgUnifiedfS1(y,dydx,x,htry,eps,yscal,hdid,hnext)
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
      double precision dydx(5,0:nint_max)
      double precision y(5,0:nint_max)
      double precision yscal(5,0:nint_max)
      double precision errmax,h,htemp,xnew
      double precision yerr(5,0:nint_max)
      double precision ytemp(5,0:nint_max)
*
      h = htry
*
 101  call rkcksgUnifiedfS1(y,dydx,x,h,ytemp,yerr)
*
      errmax = 0d0
      do i=1,5
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
         do i=1,5
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
      subroutine rkcksgUnifiedfS1(y,dydx,x,h,yout,yerr)
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
      double precision dydx(5,0:nint_max)
      double precision y(5,0:nint_max)
      double precision yerr(5,0:nint_max)
      double precision yout(5,0:nint_max)
      double precision ytemp(5,0:nint_max)
      double precision ak2(5,0:nint_max)
      double precision ak3(5,0:nint_max)
      double precision ak4(5,0:nint_max)
      double precision ak5(5,0:nint_max)
      double precision ak6(5,0:nint_max)
*
      do i=1,5
         do alpha=0,nin(igrid)
            ytemp(i,alpha) = y(i,alpha)
     1                     + B21 * h * dydx(i,alpha)
         enddo
      enddo
      call derivssgUnifiedfS1(x+A2*h,ytemp,ak2)
      do i=1,5
         do alpha=0,nin(igrid)
            ytemp(i,alpha) = y(i,alpha) 
     1                     + h * ( B31 * dydx(i,alpha)
     2                     +       B32 * ak2(i,alpha) )
         enddo
      enddo
      call derivssgUnifiedfS1(x+A3*h,ytemp,ak3)
      do i=1,5
         do alpha=0,nin(igrid)
            ytemp(i,alpha) = y(i,alpha) 
     1                     + h * ( B41 * dydx(i,alpha)
     2                     +       B42 * ak2(i,alpha) 
     3                     +       B43 * ak3(i,alpha) )
         enddo
      enddo
      call derivssgUnifiedfS1(x+A4*h,ytemp,ak4)
      do i=1,5
         do alpha=0,nin(igrid)
            ytemp(i,alpha) = y(i,alpha) 
     1                     + h * ( B51 * dydx(i,alpha)
     2                     +       B52 * ak2(i,alpha)
     3                     +       B53 * ak3(i,alpha)
     4                     +       B54 * ak4(i,alpha) )
         enddo
      enddo
      call derivssgUnifiedfS1(x+A5*h,ytemp,ak5)
      do i=1,5
         do alpha=0,nin(igrid)
            ytemp(i,alpha) = y(i,alpha)
     1                     + h * ( B61 * dydx(i,alpha)
     2                     +       B62 * ak2(i,alpha)
     3                     +       B63 * ak3(i,alpha)
     4                     +       B64 * ak4(i,alpha)
     5                     +       B65 * ak5(i,alpha) )
         enddo
      enddo
      call derivssgUnifiedfS1(x+A6*h,ytemp,ak6)
      do i=1,5
         do alpha=0,nin(igrid)
            yout(i,alpha) = y(i,alpha)
     1                    + h * ( C1 * dydx(i,alpha)
     2                    +       C3 * ak3(i,alpha)
     3                    +       C4 * ak4(i,alpha)
     4                    +       C6 * ak6(i,alpha) )
         enddo
      enddo
      do i=1,5
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
      subroutine derivssgUnifiedfS1(t,f,dfdt)
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
      double precision t
      double precision f(5,0:nint_max)
**
*     Internal Variables
*
      integer i,l
      integer alpha,delta
      double precision mu2
      double precision integralsQCD
      double precision integralsQED
      double precision coupQCD,a_QCD,muR2,bts,fbeta
      double precision coupQED,a_QED
      double precision Deltaud
      double precision integ(0:nint_max,5,5)
**
*     Output Variables
*
      double precision dfdt(5,0:nint_max)
*
*     Couplings
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
*     Deltaud = ( nf_up - nf_dw ) / nf
*
      Deltaud = 0d0
      if(wnf.eq.3.or.wnf.eq.5) Deltaud = - 1d0 / dble(wnf)
*
      do alpha=0,nin(igrid)
*     QCD
         integ(alpha,1,1) = integralsQCD(0,alpha,coupQCD,7)
         integ(alpha,1,2) = 0d0
         integ(alpha,1,3) = integralsQCD(0,alpha,coupQCD,6)
         integ(alpha,1,4) = 0d0
         integ(alpha,1,5) = 0d0
*
         integ(alpha,2,1) = 0d0
         integ(alpha,2,2) = 0d0
         integ(alpha,2,3) = 0d0
         integ(alpha,2,4) = 0d0
         integ(alpha,2,5) = 0d0
*
         integ(alpha,3,1) = integralsQCD(0,alpha,coupQCD,5)
         integ(alpha,3,2) = 0d0
         integ(alpha,3,3) = integralsQCD(0,alpha,coupQCD,4)
         integ(alpha,3,4) = 0d0
         integ(alpha,3,5) = 0d0
*
         integ(alpha,4,1) = Deltaud * integralsQCD(0,alpha,coupQCD,5)
         integ(alpha,4,2) = 0d0
         integ(alpha,4,3) = Deltaud * ( integralsQCD(0,alpha,coupQCD,4)
     1                    - integralsQCD(0,alpha,coupQCD,1) )
         integ(alpha,4,4) = integralsQCD(0,alpha,coupQCD,1)
         integ(alpha,4,5) = 0d0
*
         integ(alpha,5,1) = 0d0
         integ(alpha,5,2) = 0d0
         integ(alpha,5,3) = 0d0
         integ(alpha,5,4) = 0d0
         integ(alpha,5,5) = 0d0
*     QED
         integ(alpha,2,2) = integ(alpha,2,2)
     1        + bts * integralsQED(0,alpha,coupQED,coupQCD,3)
         integ(alpha,2,3) = integ(alpha,2,3)
     1        + bts * integralsQED(0,alpha,coupQED,coupQCD,4)
         integ(alpha,2,4) = integ(alpha,2,4)
     1        + bts * integralsQED(0,alpha,coupQED,coupQCD,5)
         integ(alpha,2,5) = integ(alpha,2,5)
     1        + bts * integralsQED(0,alpha,coupQED,coupQCD,13)
*
         integ(alpha,3,2) = integ(alpha,3,2)
     1        + bts * integralsQED(0,alpha,coupQED,coupQCD,6)
         integ(alpha,3,3) = integ(alpha,3,3)
     1        + bts * integralsQED(0,alpha,coupQED,coupQCD,7)
         integ(alpha,3,4) = integ(alpha,3,4)
     1        + bts * integralsQED(0,alpha,coupQED,coupQCD,8)
*
         integ(alpha,4,2) = integ(alpha,4,2)
     1        + bts * integralsQED(0,alpha,coupQED,coupQCD,9)
         integ(alpha,4,3) = integ(alpha,4,3)
     1        + bts * integralsQED(0,alpha,coupQED,coupQCD,10)
         integ(alpha,4,4) = integ(alpha,4,4)
     1        + bts * integralsQED(0,alpha,coupQED,coupQCD,11)
*
         integ(alpha,5,2) = integ(alpha,5,2)
     1        + bts * integralsQED(0,alpha,coupQED,coupQCD,14)
         integ(alpha,5,5) = integ(alpha,5,5)
     1        + bts * integralsQED(0,alpha,coupQED,coupQCD,12)
      enddo
*     
*     Initialization
*
      do i=1,5
         do alpha=0,nin(igrid)
            dfdt(i,alpha) = 0d0
            do l=1,5
               do delta=0,nin(igrid)-alpha
                  dfdt(i,alpha) = dfdt(i,alpha)
     1                          + integ(delta,i,l) * f(l,alpha+delta)
               enddo
            enddo
         enddo
      enddo
*
      return
      end
*
************************************************************************
*
*     Second Singolet
*
************************************************************************
      subroutine odeintsgUnifiedfS2(mu21,mu22,ystart,y)
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
      double precision ystart(2,0:nint_max)
**
*     Internal Variables
*
      integer i,nstp
      integer alpha
      double precision x1,x2
      double precision a_QCD
      double precision h,hdid,hnext,x
      double precision dydx(2,0:nint_max)
      double precision yscal(2,0:nint_max)
**
*     Output Variables
*
      double precision y(2,0:nint_max)
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
      do i=1,2
         do alpha=0,nin(igrid)
            y(i,alpha) = ystart(i,alpha)
         enddo
      enddo
*
      do nstp=1,maxstp
         call derivssgUnifiedfS2(x,y,dydx)
*
         do i=1,2
            do alpha=0,nin(igrid)
               yscal(i,alpha) = dabs(y(i,alpha)) 
     1                        + dabs(h*dydx(i,alpha)) 
     2                        + tiny
            enddo
         enddo
*
         if((x+h-x2)*(x+h-x1).gt.0d0) h = x2 - x
*
         call rkqssgUnifiedfS2(y,dydx,x,h,eps,yscal,hdid,hnext)
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
      subroutine rkqssgUnifiedfS2(y,dydx,x,htry,eps,yscal,hdid,hnext)
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
      double precision dydx(2,0:nint_max)
      double precision y(2,0:nint_max)
      double precision yscal(2,0:nint_max)
      double precision errmax,h,htemp,xnew
      double precision yerr(2,0:nint_max)
      double precision ytemp(2,0:nint_max)
*
      h = htry
*
 101  call rkcksgUnifiedfS2(y,dydx,x,h,ytemp,yerr)
*
      errmax = 0d0
      do i=1,2
         do alpha=0,nin(igrid)
            errmax = max(errmax,dabs(yerr(i,alpha)
     1                              /yscal(i,alpha)))
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
         do i=1,2
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
      subroutine rkcksgUnifiedfS2(y,dydx,x,h,yout,yerr)
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
      double precision dydx(2,0:nint_max)
      double precision y(2,0:nint_max)
      double precision yerr(2,0:nint_max)
      double precision yout(2,0:nint_max)
      double precision ytemp(2,0:nint_max)
      double precision ak2(2,0:nint_max)
      double precision ak3(2,0:nint_max)
      double precision ak4(2,0:nint_max)
      double precision ak5(2,0:nint_max)
      double precision ak6(2,0:nint_max)
*
      do i=1,2
         do alpha=0,nin(igrid)
            ytemp(i,alpha) = y(i,alpha)
     1                     + B21 * h * dydx(i,alpha)
         enddo
      enddo
      call derivssgUnifiedfS2(x+A2*h,ytemp,ak2)
      do i=1,2
         do alpha=0,nin(igrid)
            ytemp(i,alpha) = y(i,alpha) 
     1                     + h * ( B31 * dydx(i,alpha)
     2                     +       B32 * ak2(i,alpha) )
         enddo
      enddo
      call derivssgUnifiedfS2(x+A3*h,ytemp,ak3)
      do i=1,2
         do alpha=0,nin(igrid)
            ytemp(i,alpha) = y(i,alpha) 
     1                     + h * ( B41 * dydx(i,alpha)
     2                     +       B42 * ak2(i,alpha) 
     3                     +       B43 * ak3(i,alpha) )
         enddo
      enddo
      call derivssgUnifiedfS2(x+A4*h,ytemp,ak4)
      do i=1,2
         do alpha=0,nin(igrid)
            ytemp(i,alpha) = y(i,alpha) 
     1                     + h * ( B51 * dydx(i,alpha)
     2                     +       B52 * ak2(i,alpha)
     3                     +       B53 * ak3(i,alpha)
     4                     +       B54 * ak4(i,alpha) )
         enddo
      enddo
      call derivssgUnifiedfS2(x+A5*h,ytemp,ak5)
      do i=1,2
         do alpha=0,nin(igrid)
            ytemp(i,alpha) = y(i,alpha)
     1                     + h * ( B61 * dydx(i,alpha)
     2                     +       B62 * ak2(i,alpha)
     3                     +       B63 * ak3(i,alpha)
     4                     +       B64 * ak4(i,alpha)
     5                     +       B65 * ak5(i,alpha) )
         enddo
      enddo
      call derivssgUnifiedfS2(x+A6*h,ytemp,ak6)
      do i=1,2
         do alpha=0,nin(igrid)
            yout(i,alpha) = y(i,alpha)
     1                    + h * ( C1 * dydx(i,alpha)
     2                    +       C3 * ak3(i,alpha)
     3                    +       C4 * ak4(i,alpha)
     4                    +       C6 * ak6(i,alpha) )
         enddo
      enddo
      do i=1,2
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
      subroutine derivssgUnifiedfS2(t,f,dfdt)
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
      double precision t
      double precision f(2,0:nint_max)
**
*     Internal Variables
*
      integer i,l
      integer alpha,delta
      double precision mu2
      double precision integralsQCD
      double precision integralsQED
      double precision coupQCD,a_QCD,muR2,bts,fbeta
      double precision coupQED,a_QED
      double precision Deltaud
      double precision integ(0:nint_max,2,2)
**
*     Output Variables
*
      double precision dfdt(2,0:nint_max)
*
*     Couplings
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
*     Deltaud = ( nf_up - nf_dw ) / nf
*
      Deltaud = 0d0
      if(wnf.eq.3.or.wnf.eq.5) Deltaud = - 1d0 / dble(wnf)
*
      do alpha=0,nin(igrid)
*     QCD
         integ(alpha,1,1) = integralsQCD(0,alpha,coupQCD,3)
         integ(alpha,1,2) = 0d0
*
         integ(alpha,2,1) = Deltaud * ( integralsQCD(0,alpha,coupQCD,3)
     1                    - integralsQCD(0,alpha,coupQCD,2) )
         integ(alpha,2,2) = integralsQCD(0,alpha,coupQCD,2)
*     QED
         integ(alpha,1,1) = integ(alpha,1,1)
     1        + bts * integralsQED(0,alpha,coupQED,coupQCD,7)
         integ(alpha,1,2) = integ(alpha,1,2)
     1        + bts * integralsQED(0,alpha,coupQED,coupQCD,8)
*
         integ(alpha,2,1) = integ(alpha,2,1)
     1        + bts * integralsQED(0,alpha,coupQED,coupQCD,10)
         integ(alpha,2,2) = integ(alpha,2,2)
     1        + bts * integralsQED(0,alpha,coupQED,coupQCD,11)
      enddo
*
*     Initialization
*
      do i=1,2
         do alpha=0,nin(igrid)
            dfdt(i,alpha) = 0d0
            do l=1,2
               do delta=0,nin(igrid)-alpha
                  dfdt(i,alpha) = dfdt(i,alpha)
     1            + integ(delta,i,l) * f(l,alpha+delta)
               enddo
            enddo
         enddo
      enddo
*
      return
      end
