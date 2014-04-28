************************************************************************
*
*     odeintsgQED.f:
*
*     Set of routines that perform the Adaptive Stepsize Control for 
*     Runge-Kutta integration of an Ordinary Differential Equation.
*
************************************************************************
      subroutine odeintsgQED(mu21,mu22,ystart,y)
*
      implicit none
*
      include "../commons/PDFEvolution.h"
      include "../commons/grid.h"
**
*     Input Variables
*
      double precision mu21,mu22
      double precision ystart(3,3,0:nint_max,0:nint_max)
**
*     Internal Variables
*
      integer i,j,nstp
      integer alpha,beta
      integer maxstp
      double precision x1,x2
      double precision a_QED
      double precision h1,eps
      double precision h,hdid,hnext,x
      double precision dydx(3,3,0:nint_max,0:nint_max)
      double precision yscal(3,3,0:nint_max,0:nint_max)
      double precision tiny

      parameter(maxstp=1000)
      parameter(tiny=1d1)
      parameter(h1=dlog(5d-1))
      parameter(eps=1d-5)
**
*     Output Variables
*
      double precision y(3,3,0:nint_max,0:nint_max)
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
         do j=1,3
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  y(i,j,alpha,beta) = ystart(i,j,alpha,beta)
               enddo
            enddo
         enddo
      enddo
*
      do nstp=1,maxstp
         call derivssgQED(x,y,dydx)
*
         do i=1,3
            do j=1,3
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
         call rkqssgQED(y,dydx,x,h,eps,yscal,hdid,hnext)
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
      subroutine rkqssgQED(y,dydx,x,htry,eps,yscal,hdid,hnext)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Variables
*
      integer i,j
      integer alpha,beta
      double precision safety,pgrow,pshrnk,errcon
      double precision eps,hdid,hnext,htry,x
      double precision dydx(3,3,0:nint_max,0:nint_max)
      double precision y(3,3,0:nint_max,0:nint_max)
      double precision yscal(3,3,0:nint_max,0:nint_max)
      double precision errmax,h,htemp,xnew
      double precision yerr(3,3,0:nint_max,0:nint_max)
      double precision ytemp(3,3,0:nint_max,0:nint_max)

      parameter(safety=0.9d0,pgrow=-0.2d0,pshrnk=-0.25d0,errcon=1.89d-4)
*
      h = htry
*
 101  call rkcksgQED(y,dydx,x,h,ytemp,yerr)
*
      errmax = 0d0
      do i=1,3
         do j=1,3
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
         do i=1,3
            do j=1,3
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
      subroutine rkcksgQED(y,dydx,x,h,yout,yerr)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Variables
*
      integer i,j
      integer alpha,beta
      double precision h,x
      double precision dydx(3,3,0:nint_max,0:nint_max)
      double precision y(3,3,0:nint_max,0:nint_max)
      double precision yerr(3,3,0:nint_max,0:nint_max)
      double precision yout(3,3,0:nint_max,0:nint_max)
      double precision ytemp(3,3,0:nint_max,0:nint_max)
      double precision ak2(3,3,0:nint_max,0:nint_max)
      double precision ak3(3,3,0:nint_max,0:nint_max)
      double precision ak4(3,3,0:nint_max,0:nint_max)
      double precision ak5(3,3,0:nint_max,0:nint_max)
      double precision ak6(3,3,0:nint_max,0:nint_max)
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
      do i=1,3
         do j=1,3
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  ytemp(i,j,alpha,beta) = y(i,j,alpha,beta)
     1                              + B21 * h * dydx(i,j,alpha,beta)
               enddo
            enddo
         enddo
      enddo
      call derivssgQED(x+A2*h,ytemp,ak2)
      do i=1,3
         do j=1,3
            do alpha=0,nin(igrid)
               do beta=0,nin(igrid)
                  ytemp(i,j,alpha,beta) = y(i,j,alpha,beta) 
     1                              + h * ( B31 * dydx(i,j,alpha,beta)
     2                              +       B32 * ak2(i,j,alpha,beta) )
               enddo
            enddo
         enddo
      enddo
      call derivssgQED(x+A3*h,ytemp,ak3)
      do i=1,3
         do j=1,3
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
      call derivssgQED(x+A4*h,ytemp,ak4)
      do i=1,3
         do j=1,3
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
      call derivssgQED(x+A5*h,ytemp,ak5)
      do i=1,3
         do j=1,3
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
      call derivssgQED(x+A6*h,ytemp,ak6)
      do i=1,3
         do j=1,3
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
      do i=1,3
         do j=1,3
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
      subroutine derivssgQED(t,M,dMdt)
*
      implicit none
*
      include "../commons/PDFEvolution.h"
      include "../commons/grid.h"
**
*     Input Variables
*
      double precision t
      double precision M(3,3,0:nint_max,0:nint_max)
**
*     Internal Variables
*
      integer i,j,k,mapp(3,3)
      integer alpha,beta,delta
      double precision mu2
      double precision integralsQED
      double precision coup,a_QED
      double precision integ(0:nint_max,3,3)
**
*     Output Variables
*
      double precision dMdt(3,3,0:nint_max,0:nint_max)
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
            do j=1,3
               integ(alpha,i,j) = integralsQED(0,alpha,coup,mapp(i,j))
            enddo
         enddo
      enddo
*
*     Initialization
*
      do i=1,3
         do j=1,3
            do alpha=0,nin(igrid)
               do beta=alpha,nin(igrid)
                  dMdt(i,j,alpha,beta) = 0d0
                  do k=1,3
                     do delta=0,nin(igrid)-alpha
                        dMdt(i,j,alpha,beta) = dMdt(i,j,alpha,beta)
     1                  + integ(delta,i,k) * M(k,j,alpha+delta,beta)
                     enddo
                  enddo
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
*     The following routine computes the derivative of the singlet
*     evolution operator in QED.
*
************************************************************************
      subroutine DeriveSgQED(coup,dMdt)
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
      integer i,j,mapp(3,3)
      integer alpha,beta
      double precision integralsQED
      double precision integ(0:nint_max,3,3)
**
*     Output Variables
*
      double precision dMdt(3,3,0:nint_max,0:nint_max)
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
            do j=1,3
               integ(alpha,i,j) = integralsQED(0,alpha,coup,mapp(i,j))
            enddo
         enddo
      enddo
*
*     Initialization
*
      do i=1,3
         do j=1,3
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
