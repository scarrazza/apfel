************************************************************************
*
*     lambda.f:
*
*     This routins returns the value of lambdaQCD with nf flavours
*     given the reference value with nf=5 active flavours.
*
************************************************************************
      subroutine lambda(nf,lref,xlam)
*
      implicit none
*
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer nf
      double precision lref
**
*     Internal Variables
*
      double precision lref2
      double precision mc,mb,mc2,mb2
      double precision xlam4
**
*     Output Variables
*
      double precision xlam
*
      xlam  = lref
      lref2 = lref**2
*
      if(nf.eq.4)then
         mb2  = m2th(5)
         mb   = dsqrt(mb2)
         xlam = lref * (lref/mb)**(-2d0/25d0)
     1        * (dlog(mb2/lref2))**(963d0/14375d0)
         xlam4 = xlam
      elseif(nf.eq.3)then
         mc2   = m2th(4)
         mc    = dsqrt(mc2)
         mb2   = m2th(5)
         mb    = dsqrt(mb2)
         xlam4 = lref * (lref/mb)**(-2d0/25d0)
     1         * (dlog(mb2/lref2))**(963d0/14375d0)
         xlam  = xlam4 * (mc/xlam4)**(2d0/25d0)
     1         * (dlog(mc2/xlam4**2))**(107d0/2025d0)
      elseif(nf.gt.5)then
         write(*,*) "In lambda.f:"
         write(*,*) "error: nf > 5!"
      endif
*
      return
      end
*
************************************************************************
*
*     The following function computes LambdaQCD given the number of 
*     active flavours nf and the values of alphas at the reference scale
*     mu2.
*
************************************************************************
      function lambdaNF(nf,alphas,mu2)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/ipt.h"
**
*     Input Variables
*
      integer nf
      double precision alphas,mu2
**
*     Internal Variables
*
      double precision as
      double precision lambda2,step
      double precision as_lambda
      double precision zero1,zero2
      double precision eps
      parameter(eps=1d-5)
**
*     Output Variables
*
      double precision lambdaNF
*
      as = alphas / 4d0 / pi
*
      lambda2 = 0.0025
      zero1   = as - as_lambda(nf,lambda2,mu2,ipt)
*
      step    = 0.04
      do
         lambda2 = lambda2 + step
         zero2   = as - as_lambda(nf,lambda2,mu2,ipt)
*
         if(zero1*zero2.lt.0d0)then
            if(dabs(zero2).le.eps) goto 102
            lambda2 = lambda2 - step
            if(dabs(zero1).le.eps) goto 102
            step    = step / 2d0
         else
            zero1 = zero2
         endif
      enddo
*
 102  lambdaNF = dsqrt(lambda2)
*
      return
      end
