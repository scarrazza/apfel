************************************************************************
*
*     odeintsgQCD.f:
*
*     Set of routines that perform the Adaptive Stepsize Control for 
*     Runge-Kutta integration of an Ordinary Differential Equation.
*
************************************************************************
      subroutine odeintsgQCD(mu21,mu22,ystart,y)
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
         call derivssgQCD(x,y,dydx)
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
         call rkqssgQCD(y,dydx,x,h,eps,yscal,hdid,hnext)
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
      subroutine rkqssgQCD(y,dydx,x,htry,eps,yscal,hdid,hnext)
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
 101  call rkcksgQCD(y,dydx,x,h,ytemp,yerr)
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
      subroutine rkcksgQCD(y,dydx,x,h,yout,yerr)
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
      call derivssgQCD(x+A2*h,ytemp,ak2)
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
      call derivssgQCD(x+A3*h,ytemp,ak3)
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
      call derivssgQCD(x+A4*h,ytemp,ak4)
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
      call derivssgQCD(x+A5*h,ytemp,ak5)
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
      call derivssgQCD(x+A6*h,ytemp,ak6)
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
      subroutine derivssgQCD(t,M,dMdt)
*
      implicit none
*
      include "../commons/PDFEvolution.h"
      include "../commons/grid.h"
**
*     Input Variables
*
      double precision t
      double precision M(2,2,0:nint_max,0:nint_max)
**
*     Internal Variables
*
      integer i,j,k,mapp(2,2)
      integer alpha,beta,delta
      double precision mu2
      double precision integralsQCD
      double precision coup,a_QCD
      double precision integ1(0:nint_max,2,2)
      double precision integ2(0:nint_max,0:nint_max,2,2)
**
*     Output Variables
*
      double precision dMdt(2,2,0:nint_max,0:nint_max)
*
      if(PDFEvol.eq."exactmu")then
         mu2  = dexp(t)
         coup = a_QCD(mu2)
      else
         coup = t
      endif
*
*     Map used for the muliplication
*
      mapp(1,1) = 4
      mapp(1,2) = 5
      mapp(2,1) = 6
      mapp(2,2) = 7
*
*     Perform the convolution between the splitting matrix and the 
*     evolution operator.
*     If the computation evolution operator has been enabled, perform
*     the convolution in a more general way (no possibility to use the
*     shift symmetry of the splitting matrix).
*
      if(IsExt(igrid))then
         do alpha=0,nin(igrid)
            do beta=alpha,nin(igrid)
               do i=1,2
                  do j=1,2
                     integ2(alpha,beta,i,j) = 
     1               integralsQCD(alpha,beta,coup,mapp(i,j))
                  enddo
               enddo
            enddo
         enddo
*
*     Convolution
*
         do i=1,2
            do j=1,2
               do alpha=0,nin(igrid)
                  do beta=alpha,nin(igrid)
                     dMdt(i,j,alpha,beta) = 0d0
                     do k=1,2
                        do delta=0,nin(igrid)
                           dMdt(i,j,alpha,beta) = dMdt(i,j,alpha,beta)
     1                     + integ2(alpha,delta,i,k) * M(k,j,delta,beta)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else
         do alpha=0,nin(igrid)
            do i=1,2
               do j=1,2
                  integ1(alpha,i,j) = 
     1            integralsQCD(0,alpha,coup,mapp(i,j))
               enddo
            enddo
         enddo
*
*     Convolution
*
         do i=1,2
            do j=1,2
               do alpha=0,nin(igrid)
                  do beta=alpha,nin(igrid)
                     dMdt(i,j,alpha,beta) = 0d0
                     do k=1,2
                        do delta=0,nin(igrid)-alpha
                           dMdt(i,j,alpha,beta) = dMdt(i,j,alpha,beta)
     1                     + integ1(delta,i,k) * M(k,j,alpha+delta,beta)
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
*     The following routine computes the derivative of the singlet
*     evolution operator in QCD.
*
************************************************************************
      subroutine DeriveSgQCD(coup,dMdt)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      double precision coup
**
*     Internal Variables
*
      integer i,j,mapp(2,2)
      integer alpha,beta
      double precision integralsQCD
      double precision integ(0:nint_max,2,2)
**
*     Output Variables
*
      double precision dMdt(2,2,0:nint_max,0:nint_max)
*
*     Map used for the muliplication
*
      mapp(1,1) = 4
      mapp(1,2) = 5
      mapp(2,1) = 6
      mapp(2,2) = 7
*
      do alpha=0,nin(igrid)
         do i=1,2
            do j=1,2
               integ(alpha,i,j) = integralsQCD(0,alpha,coup,mapp(i,j))
            enddo
         enddo
      enddo
*
*     Initialization
*
      do i=1,2
         do j=1,2
            do alpha=0,nin(igrid)
               do beta=0,alpha-1
                  dMdt(i,j,alpha,beta) = 0d0
               enddo
               do beta=alpha,nin(igrid)
                  dMdt(i,j,alpha,beta) = integ(beta-alpha,i,j)
               enddo
            enddo
         enddo
      enddo
*
      return
      end
