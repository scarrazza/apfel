************************************************************************
*
*     odeintsgUnified.f:
*
*     Set of routines that perform the Adaptive Stepsize Control for 
*     Runge-Kutta integration of an Ordinary Differential Equation.
*
************************************************************************
*
*     First Singolet
*
************************************************************************
      subroutine odeintsgUnifiedS1(mu21,mu22,ystart,y)
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
      double precision ystart(5,5,0:nint_max,0:nint_max)
**
*     Internal Variables
*
      integer i,j,nstp
      integer alpha,beta
      double precision x1,x2
      double precision a_QCD
      double precision h,hdid,hnext,x
      double precision dydx(5,5,0:nint_max,0:nint_max)
      double precision yscal(5,5,0:nint_max,0:nint_max)
**
*     Output Variables
*
      double precision y(5,5,0:nint_max,0:nint_max)
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
         do j=1,5
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  y(i,j,alpha,beta) = ystart(i,j,alpha,beta)
               enddo
            enddo
         enddo
      enddo
*
      do nstp=1,maxstp
         call derivssgUnifiedS1(x,y,dydx)
*
         do i=1,5
            do j=1,5
               do alpha=0,nin(igrid)
                  do beta=0,nin(igrid)
                     yscal(i,j,alpha,beta) = dabs(y(i,j,alpha,beta)) 
     1                                    + dabs(h*dydx(i,j,alpha,beta)) 
     2                                    + tiny
                  enddo
               enddo
            enddo
         enddo
*
         if((x+h-x2)*(x+h-x1).gt.0d0) h = x2 - x
*
         call rkqssgUnifiedS1(y,dydx,x,h,eps,yscal,hdid,hnext)
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
      subroutine rkqssgUnifiedS1(y,dydx,x,htry,eps,yscal,hdid,hnext)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/odeint2.h"
**
*     Variables
*
      integer i,j
      integer alpha,beta
      double precision eps,hdid,hnext,htry,x
      double precision dydx(5,5,0:nint_max,0:nint_max)
      double precision y(5,5,0:nint_max,0:nint_max)
      double precision yscal(5,5,0:nint_max,0:nint_max)
      double precision errmax,h,htemp,xnew
      double precision yerr(5,5,0:nint_max,0:nint_max)
      double precision ytemp(5,5,0:nint_max,0:nint_max)
*
      h = htry
*
 101  call rkcksgUnifiedS1(y,dydx,x,h,ytemp,yerr)
*
      errmax = 0d0
      do i=1,5
         do j=1,5
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  errmax = max(errmax,dabs(yerr(i,j,alpha,beta)
     1                                    /yscal(i,j,alpha,beta)))
               enddo
            enddo
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
            do j=1,5
               do alpha=0,nin(igrid)
                  do beta=0,nin(igrid)
                     y(i,j,alpha,beta) = ytemp(i,j,alpha,beta)
                  enddo
               enddo
            enddo
         enddo
         return
      endif
*
      end
*
************************************************************************
      subroutine rkcksgUnifiedS1(y,dydx,x,h,yout,yerr)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/odeint2.h"
**
*     Variables
*
      integer i,j
      integer alpha,beta
      double precision h,x
      double precision dydx(5,5,0:nint_max,0:nint_max)
      double precision y(5,5,0:nint_max,0:nint_max)
      double precision yerr(5,5,0:nint_max,0:nint_max)
      double precision yout(5,5,0:nint_max,0:nint_max)
      double precision ytemp(5,5,0:nint_max,0:nint_max)
      double precision ak2(5,5,0:nint_max,0:nint_max)
      double precision ak3(5,5,0:nint_max,0:nint_max)
      double precision ak4(5,5,0:nint_max,0:nint_max)
      double precision ak5(5,5,0:nint_max,0:nint_max)
      double precision ak6(5,5,0:nint_max,0:nint_max)
*
      do i=1,5
         do j=1,5
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  ytemp(i,j,alpha,beta) = y(i,j,alpha,beta)
     1                              + B21 * h * dydx(i,j,alpha,beta)
               enddo
            enddo
         enddo
      enddo
      call derivssgUnifiedS1(x+A2*h,ytemp,ak2)
      do i=1,5
         do j=1,5
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  ytemp(i,j,alpha,beta) = y(i,j,alpha,beta) 
     1                              + h * ( B31 * dydx(i,j,alpha,beta)
     2                              +       B32 * ak2(i,j,alpha,beta) )
               enddo
            enddo
         enddo
      enddo
      call derivssgUnifiedS1(x+A3*h,ytemp,ak3)
      do i=1,5
         do j=1,5
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  ytemp(i,j,alpha,beta) = y(i,j,alpha,beta) 
     1                               + h * ( B41 * dydx(i,j,alpha,beta)
     2                               +       B42 * ak2(i,j,alpha,beta) 
     3                               +       B43 * ak3(i,j,alpha,beta) )
               enddo
            enddo
         enddo
      enddo
      call derivssgUnifiedS1(x+A4*h,ytemp,ak4)
      do i=1,5
         do j=1,5
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  ytemp(i,j,alpha,beta) = y(i,j,alpha,beta) 
     1                               + h * ( B51 * dydx(i,j,alpha,beta)
     2                               +       B52 * ak2(i,j,alpha,beta)
     3                               +       B53 * ak3(i,j,alpha,beta)
     4                               +       B54 * ak4(i,j,alpha,beta) )
               enddo
            enddo
         enddo
      enddo
      call derivssgUnifiedS1(x+A5*h,ytemp,ak5)
      do i=1,5
         do j=1,5
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  ytemp(i,j,alpha,beta) = y(i,j,alpha,beta)
     1                               + h * ( B61 * dydx(i,j,alpha,beta)
     2                               +       B62 * ak2(i,j,alpha,beta)
     3                               +       B63 * ak3(i,j,alpha,beta)
     4                               +       B64 * ak4(i,j,alpha,beta)
     5                               +       B65 * ak5(i,j,alpha,beta) )
               enddo
            enddo
         enddo
      enddo
      call derivssgUnifiedS1(x+A6*h,ytemp,ak6)
      do i=1,5
         do j=1,5
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  yout(i,j,alpha,beta) = y(i,j,alpha,beta)
     1                               + h * ( C1 * dydx(i,j,alpha,beta)
     2                               +       C3 * ak3(i,j,alpha,beta)
     3                               +       C4 * ak4(i,j,alpha,beta)
     4                               +       C6 * ak6(i,j,alpha,beta) )
               enddo
            enddo
         enddo
      enddo
      do i=1,5
         do j=1,5
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  yerr(i,j,alpha,beta) = h *( DC1 * dydx(i,j,alpha,beta)
     1                               +      DC3 * ak3(i,j,alpha,beta)
     2                               +      DC4 * ak4(i,j,alpha,beta)
     3                               +      DC5 * ak5(i,j,alpha,beta)
     4                               +      DC6 * ak6(i,j,alpha,beta) )
               enddo
            enddo
         enddo
      enddo
*
      return
      end
*
************************************************************************
      subroutine derivssgUnifiedS1(t,M,dMdt)
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
      double precision M(5,5,0:nint_max,0:nint_max)
**
*     Internal Variables
*
      integer i,j,l
      integer alpha,beta,delta
      double precision mu2
      double precision integralsQCD
      double precision integralsQED
      double precision coupQCD,a_QCD,muR2,bts,fbeta
      double precision coupQED,a_QED
      double precision Deltaud
      double precision integ1(0:nint_max,5,5)
      double precision integ2(0:nint_max,0:nint_max,5,5)
**
*     Output Variables
*
      double precision dMdt(5,5,0:nint_max,0:nint_max)
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
      if(IsExt(igrid))then
         do alpha=0,nin(igrid)
            do beta=alpha,nin(igrid)
*     QCD
            integ2(alpha,beta,1,1) = integralsQCD(alpha,beta,coupQCD,7)
            integ2(alpha,beta,1,2) = 0d0
            integ2(alpha,beta,1,3) = integralsQCD(alpha,beta,coupQCD,6)
            integ2(alpha,beta,1,4) = 0d0
            integ2(alpha,beta,1,5) = 0d0
*
            integ2(alpha,beta,2,1) = 0d0
            integ2(alpha,beta,2,2) = 0d0
            integ2(alpha,beta,2,3) = 0d0
            integ2(alpha,beta,2,4) = 0d0
            integ2(alpha,beta,2,5) = 0d0
*
            integ2(alpha,beta,3,1) = integralsQCD(alpha,beta,coupQCD,5)
            integ2(alpha,beta,3,2) = 0d0
            integ2(alpha,beta,3,3) = integralsQCD(alpha,beta,coupQCD,4)
            integ2(alpha,beta,3,4) = 0d0
            integ2(alpha,beta,3,5) = 0d0
*
            integ2(alpha,beta,4,1) = Deltaud
     1           * integralsQCD(alpha,beta,coupQCD,5)
            integ2(alpha,beta,4,2) = 0d0
            integ2(alpha,beta,4,3) = Deltaud
     1           * ( integralsQCD(alpha,beta,coupQCD,4)
     2           - integralsQCD(alpha,beta,coupQCD,1) )
            integ2(alpha,beta,4,4) = integralsQCD(alpha,beta,coupQCD,1)
            integ2(alpha,beta,4,5) = 0d0
*
            integ2(alpha,beta,5,1) = 0d0
            integ2(alpha,beta,5,2) = 0d0
            integ2(alpha,beta,5,3) = 0d0
            integ2(alpha,beta,5,4) = 0d0
            integ2(alpha,beta,5,5) = 0d0
*     QED
            integ2(alpha,beta,1,1) = integ2(alpha,beta,1,1)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,5)
            integ2(alpha,beta,1,2) = integ2(alpha,beta,1,2)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,6)
            integ2(alpha,beta,1,3) = integ2(alpha,beta,1,3)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,7)
            integ2(alpha,beta,1,4) = integ2(alpha,beta,1,4)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,8)
*
            integ2(alpha,beta,2,1) = integ2(alpha,beta,2,1)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,9)
            integ2(alpha,beta,2,2) = integ2(alpha,beta,2,2)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,10)
            integ2(alpha,beta,2,3) = integ2(alpha,beta,2,3)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,11)
            integ2(alpha,beta,2,4) = integ2(alpha,beta,2,4)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,12)
            integ2(alpha,beta,2,5) = integ2(alpha,beta,2,5)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,24)
*
            integ2(alpha,beta,3,1) = integ2(alpha,beta,3,1)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,13)
            integ2(alpha,beta,3,2) = integ2(alpha,beta,3,2)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,14)
            integ2(alpha,beta,3,3) = integ2(alpha,beta,3,3)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,15)
            integ2(alpha,beta,3,4) = integ2(alpha,beta,3,4)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,16)
*
            integ2(alpha,beta,4,1) = integ2(alpha,beta,4,1)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,17)
            integ2(alpha,beta,4,2) = integ2(alpha,beta,4,2)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,18)
            integ2(alpha,beta,4,3) = integ2(alpha,beta,4,3)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,19)
            integ2(alpha,beta,4,4) = integ2(alpha,beta,4,4)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,20)
*
            integ2(alpha,beta,5,2) = integ2(alpha,beta,5,2)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,25)
            integ2(alpha,beta,5,5) = integ2(alpha,beta,5,5)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,23)
            enddo
         enddo
*
*     Initialization
*
         do i=1,5
            do j=1,5
               do alpha=0,nin(igrid)
                  do beta=alpha,nin(igrid)
                     dMdt(i,j,alpha,beta) = 0d0
                     do l=1,5
                        do delta=alpha,beta
                           dMdt(i,j,alpha,beta) = dMdt(i,j,alpha,beta)
     1                          + integ2(alpha,delta,i,l)
     2                          * M(l,j,delta,beta)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else
         do alpha=0,nin(igrid)
*     QCD
            integ1(alpha,1,1) = integralsQCD(0,alpha,coupQCD,7)
            integ1(alpha,1,2) = 0d0
            integ1(alpha,1,3) = integralsQCD(0,alpha,coupQCD,6)
            integ1(alpha,1,4) = 0d0
            integ1(alpha,1,5) = 0d0
*
            integ1(alpha,2,1) = 0d0
            integ1(alpha,2,2) = 0d0
            integ1(alpha,2,3) = 0d0
            integ1(alpha,2,4) = 0d0
            integ1(alpha,2,5) = 0d0
*
            integ1(alpha,3,1) = integralsQCD(0,alpha,coupQCD,5)
            integ1(alpha,3,2) = 0d0
            integ1(alpha,3,3) = integralsQCD(0,alpha,coupQCD,4)
            integ1(alpha,3,4) = 0d0
            integ1(alpha,3,5) = 0d0
*
            integ1(alpha,4,1) = Deltaud
     1           * integralsQCD(0,alpha,coupQCD,5)
            integ1(alpha,4,2) = 0d0
            integ1(alpha,4,3) = Deltaud
     1           * ( integralsQCD(0,alpha,coupQCD,4)
     2           - integralsQCD(0,alpha,coupQCD,1) )
            integ1(alpha,4,4) = integralsQCD(0,alpha,coupQCD,1)
            integ1(alpha,4,5) = 0d0
*
            integ1(alpha,5,1) = 0d0
            integ1(alpha,5,2) = 0d0
            integ1(alpha,5,3) = 0d0
            integ1(alpha,5,4) = 0d0
            integ1(alpha,5,5) = 0d0
*     QED
            integ1(alpha,1,1) = integ1(alpha,1,1)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,5)
            integ1(alpha,1,2) = integ1(alpha,1,2)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,6)
            integ1(alpha,1,3) = integ1(alpha,1,3)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,7)
            integ1(alpha,1,4) = integ1(alpha,1,4)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,8)
*
            integ1(alpha,2,1) = integ1(alpha,2,1)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,9)
            integ1(alpha,2,2) = integ1(alpha,2,2)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,10)
            integ1(alpha,2,3) = integ1(alpha,2,3)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,11)
            integ1(alpha,2,4) = integ1(alpha,2,4)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,12)
            integ1(alpha,2,5) = integ1(alpha,2,5)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,24)
*
            integ1(alpha,3,1) = integ1(alpha,3,1)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,13)
            integ1(alpha,3,2) = integ1(alpha,3,2)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,14)
            integ1(alpha,3,3) = integ1(alpha,3,3)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,15)
            integ1(alpha,3,4) = integ1(alpha,3,4)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,16)
*
            integ1(alpha,4,1) = integ1(alpha,4,1)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,17)
            integ1(alpha,4,2) = integ1(alpha,4,2)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,18)
            integ1(alpha,4,3) = integ1(alpha,4,3)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,19)
            integ1(alpha,4,4) = integ1(alpha,4,4)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,20)
*
            integ1(alpha,5,2) = integ1(alpha,5,2)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,25)
            integ1(alpha,5,5) = integ1(alpha,5,5)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,23)
         enddo
*
*     Initialization
*
         do i=1,5
            do j=1,5
               do alpha=0,nin(igrid)
                  do beta=alpha,nin(igrid)
                     dMdt(i,j,alpha,beta) = 0d0
                     do l=1,5
                        do delta=0,nin(igrid)-alpha
                           dMdt(i,j,alpha,beta) = dMdt(i,j,alpha,beta)
     1                          + integ1(delta,i,l)
     2                          * M(l,j,alpha+delta,beta)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
*
      return
      end
*
************************************************************************
*
*     Second Singolet
*
************************************************************************
      subroutine odeintsgUnifiedS2(mu21,mu22,ystart,y)
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
      double precision ystart(2,2,0:nint_max,0:nint_max)
**
*     Internal Variables
*
      integer i,j,nstp
      integer alpha,beta
      double precision x1,x2
      double precision a_QCD
      double precision h,hdid,hnext,x
      double precision dydx(2,2,0:nint_max,0:nint_max)
      double precision yscal(2,2,0:nint_max,0:nint_max)
**
*     Output Variables
*
      double precision y(2,2,0:nint_max,0:nint_max)
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
         do j=1,2
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  y(i,j,alpha,beta) = ystart(i,j,alpha,beta)
               enddo
            enddo
         enddo
      enddo
*
      do nstp=1,maxstp
         call derivssgUnifiedS2(x,y,dydx)
*
         do i=1,2
            do j=1,2
               do alpha=0,nin(igrid)
                  do beta=0,nin(igrid)
                     yscal(i,j,alpha,beta) = dabs(y(i,j,alpha,beta)) 
     1                                    + dabs(h*dydx(i,j,alpha,beta)) 
     2                                    + tiny
                  enddo
               enddo
            enddo
         enddo
*
         if((x+h-x2)*(x+h-x1).gt.0d0) h = x2 - x
*
         call rkqssgUnifiedS2(y,dydx,x,h,eps,yscal,hdid,hnext)
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
      subroutine rkqssgUnifiedS2(y,dydx,x,htry,eps,yscal,hdid,hnext)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/odeint2.h"
**
*     Variables
*
      integer i,j
      integer alpha,beta
      double precision eps,hdid,hnext,htry,x
      double precision dydx(2,2,0:nint_max,0:nint_max)
      double precision y(2,2,0:nint_max,0:nint_max)
      double precision yscal(2,2,0:nint_max,0:nint_max)
      double precision errmax,h,htemp,xnew
      double precision yerr(2,2,0:nint_max,0:nint_max)
      double precision ytemp(2,2,0:nint_max,0:nint_max)
*
      h = htry
*
 101  call rkcksgUnifiedS2(y,dydx,x,h,ytemp,yerr)
*
      errmax = 0d0
      do i=1,2
         do j=1,2
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  errmax = max(errmax,dabs(yerr(i,j,alpha,beta)
     1                                    /yscal(i,j,alpha,beta)))
               enddo
            enddo
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
            do j=1,2
               do alpha=0,nin(igrid)
                  do beta=0,nin(igrid)
                     y(i,j,alpha,beta) = ytemp(i,j,alpha,beta)
                  enddo
               enddo
            enddo
         enddo
         return
      endif
*
      end
*
************************************************************************
      subroutine rkcksgUnifiedS2(y,dydx,x,h,yout,yerr)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/odeint2.h"
**
*     Variables
*
      integer i,j
      integer alpha,beta
      double precision h,x
      double precision dydx(2,2,0:nint_max,0:nint_max)
      double precision y(2,2,0:nint_max,0:nint_max)
      double precision yerr(2,2,0:nint_max,0:nint_max)
      double precision yout(2,2,0:nint_max,0:nint_max)
      double precision ytemp(2,2,0:nint_max,0:nint_max)
      double precision ak2(2,2,0:nint_max,0:nint_max)
      double precision ak3(2,2,0:nint_max,0:nint_max)
      double precision ak4(2,2,0:nint_max,0:nint_max)
      double precision ak5(2,2,0:nint_max,0:nint_max)
      double precision ak6(2,2,0:nint_max,0:nint_max)
*
      do i=1,2
         do j=1,2
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  ytemp(i,j,alpha,beta) = y(i,j,alpha,beta)
     1                              + B21 * h * dydx(i,j,alpha,beta)
               enddo
            enddo
         enddo
      enddo
      call derivssgUnifiedS2(x+A2*h,ytemp,ak2)
      do i=1,2
         do j=1,2
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  ytemp(i,j,alpha,beta) = y(i,j,alpha,beta) 
     1                              + h * ( B31 * dydx(i,j,alpha,beta)
     2                              +       B32 * ak2(i,j,alpha,beta) )
               enddo
            enddo
         enddo
      enddo
      call derivssgUnifiedS2(x+A3*h,ytemp,ak3)
      do i=1,2
         do j=1,2
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  ytemp(i,j,alpha,beta) = y(i,j,alpha,beta) 
     1                               + h * ( B41 * dydx(i,j,alpha,beta)
     2                               +       B42 * ak2(i,j,alpha,beta) 
     3                               +       B43 * ak3(i,j,alpha,beta) )
               enddo
            enddo
         enddo
      enddo
      call derivssgUnifiedS2(x+A4*h,ytemp,ak4)
      do i=1,2
         do j=1,2
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  ytemp(i,j,alpha,beta) = y(i,j,alpha,beta) 
     1                               + h * ( B51 * dydx(i,j,alpha,beta)
     2                               +       B52 * ak2(i,j,alpha,beta)
     3                               +       B53 * ak3(i,j,alpha,beta)
     4                               +       B54 * ak4(i,j,alpha,beta) )
               enddo
            enddo
         enddo
      enddo
      call derivssgUnifiedS2(x+A5*h,ytemp,ak5)
      do i=1,2
         do j=1,2
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  ytemp(i,j,alpha,beta) = y(i,j,alpha,beta)
     1                               + h * ( B61 * dydx(i,j,alpha,beta)
     2                               +       B62 * ak2(i,j,alpha,beta)
     3                               +       B63 * ak3(i,j,alpha,beta)
     4                               +       B64 * ak4(i,j,alpha,beta)
     5                               +       B65 * ak5(i,j,alpha,beta) )
               enddo
            enddo
         enddo
      enddo
      call derivssgUnifiedS2(x+A6*h,ytemp,ak6)
      do i=1,2
         do j=1,2
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  yout(i,j,alpha,beta) = y(i,j,alpha,beta)
     1                               + h * ( C1 * dydx(i,j,alpha,beta)
     2                               +       C3 * ak3(i,j,alpha,beta)
     3                               +       C4 * ak4(i,j,alpha,beta)
     4                               +       C6 * ak6(i,j,alpha,beta) )
               enddo
            enddo
         enddo
      enddo
      do i=1,2
         do j=1,2
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  yerr(i,j,alpha,beta) = h *( DC1 * dydx(i,j,alpha,beta)
     1                               +      DC3 * ak3(i,j,alpha,beta)
     2                               +      DC4 * ak4(i,j,alpha,beta)
     3                               +      DC5 * ak5(i,j,alpha,beta)
     4                               +      DC6 * ak6(i,j,alpha,beta) )
               enddo
            enddo
         enddo
      enddo
*
      return
      end
*
************************************************************************
      subroutine derivssgUnifiedS2(t,M,dMdt)
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
      double precision M(2,2,0:nint_max,0:nint_max)
**
*     Internal Variables
*
      integer i,j,l
      integer alpha,beta,delta
      double precision mu2
      double precision integralsQCD
      double precision integralsQED
      double precision coupQCD,a_QCD,muR2,bts,fbeta
      double precision coupQED,a_QED
      double precision Deltaud
      double precision integ1(0:nint_max,2,2)
      double precision integ2(0:nint_max,0:nint_max,2,2)
**
*     Output Variables
*
      double precision dMdt(2,2,0:nint_max,0:nint_max)
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
      if(IsExt(igrid))then
         do alpha=0,nin(igrid)
            do beta=alpha,nin(igrid)
*     QCD
            integ2(alpha,beta,1,1) = integralsQCD(alpha,beta,coupQCD,3)
            integ2(alpha,beta,1,2) = 0d0
*
            integ2(alpha,beta,2,1) = Deltaud
     1           * ( integralsQCD(alpha,beta,coupQCD,3)
     2           - integralsQCD(alpha,beta,coupQCD,2) )
            integ2(alpha,beta,2,2) = integralsQCD(alpha,beta,coupQCD,2)
*     QED
            integ2(alpha,beta,1,1) = integ2(alpha,beta,1,1)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,21)
            integ2(alpha,beta,1,2) = integ2(alpha,beta,1,2)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,22)
*
            integ2(alpha,beta,2,1) = integ2(alpha,beta,2,1)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,22)
            integ2(alpha,beta,2,2) = integ2(alpha,beta,2,2)
     1           + bts * integralsQED(alpha,beta,coupQED,coupQCD,21)
            enddo
         enddo
*
*     Initialization
*
         do i=1,2
            do j=1,2
               do alpha=0,nin(igrid)
                  do beta=alpha,nin(igrid)
                     dMdt(i,j,alpha,beta) = 0d0
                     do l=1,2
                        do delta=alpha,beta
                           dMdt(i,j,alpha,beta) = dMdt(i,j,alpha,beta)
     1                          + integ2(alpha,delta,i,l)
     2                          * M(l,j,delta,beta)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else
         do alpha=0,nin(igrid)
*     QCD
            integ1(alpha,1,1) = integralsQCD(0,alpha,coupQCD,3)
            integ1(alpha,1,2) = 0d0
*
            integ1(alpha,2,1) = Deltaud
     1           * ( integralsQCD(0,alpha,coupQCD,3)
     2           - integralsQCD(0,alpha,coupQCD,2) )
            integ1(alpha,2,2) = integralsQCD(0,alpha,coupQCD,2)
*     QED
            integ1(alpha,1,1) = integ1(alpha,1,1)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,21)
            integ1(alpha,1,2) = integ1(alpha,1,2)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,22)
*
            integ1(alpha,2,1) = integ1(alpha,2,1)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,22)
            integ1(alpha,2,2) = integ1(alpha,2,2)
     1           + bts * integralsQED(0,alpha,coupQED,coupQCD,21)
         enddo
*
*     Initialization
*
         do i=1,2
            do j=1,2
               do alpha=0,nin(igrid)
                  do beta=alpha,nin(igrid)
                     dMdt(i,j,alpha,beta) = 0d0
                     do l=1,2
                        do delta=0,nin(igrid)-alpha
                           dMdt(i,j,alpha,beta) = dMdt(i,j,alpha,beta)
     1                          + integ1(delta,i,l)
     2                          * M(l,j,alpha+delta,beta)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
*
      return
      end
