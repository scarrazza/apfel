************************************************************************
*
*     a_QCD.f:
*
*     It returns a_QCD = alpha_QCD / 4 / pi.
*
*     Note that the only argument of this routine is the factorization 
*     scale mu2F which is convereted into the renormalization scale 
*     inside the routine itself.
*
************************************************************************
      function a_QCD(mu2F)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/alpha_ref_QCD.h"
      include "../commons/lambda_ref_QCD.h"
      include "../commons/AlphaEvolution.h"
      include "../commons/kren.h"
      include "../commons/m2th.h"
      include "../commons/mass_scheme.h"
      include "../commons/Nf_FF.h"
      include "../commons/ipt.h"
      include "../commons/Evs.h"
      include "../commons/MaxFlavourAlpha.h"
**
*     Input Variables
*
      double precision mu2F
**
*     Internal Variables
*
      integer i
      integer nfi,nff
      integer dnf,snf
      double precision mur2th(4:6),lambda2(3:6)
      double precision mur2,mur20
      double precision asi,asr0
      double precision c1,c2,kappa,ln
      double precision as_exact,as_expanded,as_lambda
      external as_exact,as_expanded,as_lambda
**
*     Output Variables
*
      double precision a_QCD
*
c      a_QCD = 0d0
c      return
      asr0 = alpha_ref_QCD / 4d0 / pi
*
      mur20 = q2_ref_QCD
      mur2  = kren * mu2F
      do i=4,6
         mur2th(i) = kren * m2th(i)
c         mur2th(i) = m2th(i)
      enddo
      do i=3,6
         lambda2(i) = LambdaQCD(i)**2
      enddo
*
      if(Evs.eq."FF")then
         nfi = Nf_FF
         nff = Nf_FF
      elseif(Evs.eq."VF")then
         if(mur2.ge.mur2th(6))then
            nff = 6
         elseif(mur2.ge.mur2th(5))then
            nff = 5
         elseif(mur2.ge.mur2th(4))then
            nff = 4
         else
            nff = 3
         endif
         if(nff.gt.nfMaxAlpha) nff = nfMaxAlpha
*
         if(mur20.gt.mur2th(6))then
            nfi = 6
         elseif(mur20.gt.mur2th(5))then
            nfi = 5
         elseif(mur20.gt.mur2th(4))then
            nfi = 4
         else
            nfi = 3
         endif
         if(nfi.gt.nfMaxAlpha) nfi = nfMaxAlpha
      endif
*
 10   if(nff.eq.nfi) then
         if(AlphaEvol(1:5).eq."exact")
     1        a_QCD = as_exact(nfi,mur20,asr0,mur2,ipt)
         if(AlphaEvol(1:8).eq."expanded")
     1        a_QCD = as_expanded(nfi,mur20,asr0,mur2,ipt)
         if(AlphaEvol(1:6).eq."lambda")
     1        a_QCD = as_lambda(nfi,lambda2(nfi),mur2,ipt)
         return
      else
         if(nff.gt.nfi)then
            dnf = 1
            snf = 1
         else
            dnf = -1
            snf = 0
         endif
*
         if(AlphaEvol(1:5).eq."exact")
     1        asi = as_exact(nfi,mur20,asr0,mur2th(nfi+snf),ipt)
         if(AlphaEvol(1:8).eq."expanded")
     1        asi = as_expanded(nfi,mur20,asr0,mur2th(nfi+snf),ipt)
         if(AlphaEvol(1:6).eq."lambda")
     1        asi = as_lambda(nfi,lambda2(nfi),mur2th(nfi+snf),ipt)
*
*     c1 and c2 are the same coefficients used in eq. (2.42) of hep-ph/0408244 and 
*     obtained in eq. (10) of hep-ph/9706430. In the following they are divided by 
*     (4*pi) and (4*pi)^2 respectively to match the notations. Note that in terms 
*     of the MSbar mass this coefficients change.
*
         kappa = kren * k2th(nfi+snf)          ! mu_R / mu_F
c         kappa = kren                         ! mu_R / mu_F
c         kappa = 1d0                          ! mu_R / mu_F
         ln = dlog(kappa)
*     Pole Mass
         if(mass_scheme.eq."Pole")then
            if(nff.gt.nfi)then
               c1 = 2d0 / 3d0 * ln
               c2 = 4d0 / 9d0 * ln**2 + 38d0 / 3d0 * ln + 14d0 / 3d0
            elseif(nff.lt.nfi)then
               c1 = - 2d0 / 3d0 * ln
               c2 = 4d0 / 9d0 * ln**2 - 38d0 / 3d0 * ln - 14d0 / 3d0
            endif
*     MSbar mass
         elseif(mass_scheme.eq."MSbar")then
            if(nff.gt.nfi)then
               c1 = 2d0 / 3d0 * ln
               c2 = 4d0 / 9d0 * ln**2 + 22d0 / 3d0 * ln - 22d0 / 9d0
            elseif(nff.lt.nfi)then
               c1 = - 2d0 / 3d0 * ln
               c2 = 4d0 / 9d0 * ln**2 - 22d0 / 3d0 * ln + 22d0 / 9d0
            endif
         endif
*
*     NLO and NNLO threshold matchings
*
         if(ipt.eq.1)then
            asi = asi * ( 1d0 + c1 * asi )
         elseif(ipt.eq.2)then
            asi = asi * ( 1d0 + c1 * asi + c2 * asi**2 )
         endif
*
         asr0  = asi
         mur20 = mur2th(nfi+snf)
         nfi  = nfi + dnf
         goto 10
      endif
*
      end
*
************************************************************************
*
*     Routines for the computation of alpha_s with fixed number of 
*     flavours
* 
*     - as_expanded: computes alpha_s as function of alpha_s at a given
*                    refernce scale.
*
*     - as_exact: exact solution of the QCD beta function equation using
*                 fourth order Runge-Kutta algorithm.
*
*     - as_lambda: solution in terms of LambdaQCD.
*
************************************************************************
      function as_expanded(nf,mu20,as0,mu2,ipt)
*
      implicit none
**
*     Input Variables
*
      integer nf,ipt
      double precision mu2
      double precision mu20,as0
**
*     Internal Variables
*
      double precision asi
      double precision alo,t,as,den
      double precision beta0apf,beta1apf,beta2apf,b1,b2
**
*     Output Variables
*
      double precision as_expanded
*
      b1 = beta1apf(nf) / beta0apf(nf)
      b2 = beta2apf(nf) / beta0apf(nf)
*
      asi = as0
      t   = log(mu2/mu20)
      den = 1d0 + beta0apf(nf) * asi * t
      alo = asi / den
*
*     LO
*
      as = alo
*
*     NLO
*
      if(ipt.ge.1)then
         as = alo * ( 1d0 - b1 * alo * log(den) )
      endif
*
*     NNLO
*
      if(ipt.eq.2)then
         as = alo * ( 1d0 
     1      + ( alo * ( alo - asi ) * ( b2 - b1**2 )
     2      + as * b1 * dlog(as/asi) ) )
      endif
*
      as_expanded = as
*
      return
      end
*
****************************************************************************
      FUNCTION AS_EXACT(NF,MU20,AS0,MU2,IPT)
*
      IMPLICIT NONE
**
*     Input Variables
*
      INTEGER NF,IPT
      DOUBLE PRECISION MU2
      DOUBLE PRECISION AS0,MU20
**
*     Internal Variables
*
      INTEGER NSTEP,K1
      DOUBLE PRECISION AS
      DOUBLE PRECISION BETA0APF,FBETA
      DOUBLE PRECISION DLR,LRRAT,SXTH
      DOUBLE PRECISION XK0,XK1,XK2,XK3
      PARAMETER(NSTEP=10)
      PARAMETER(SXTH=0.166666666666666D0)
**
*     Output Variables
*
      DOUBLE PRECISION AS_EXACT
*
      AS    = AS0
      LRRAT = DLOG( MU2 / MU20 )
      DLR   = LRRAT / NSTEP
*
*     Analytical solution at leading order
*
      IF(IPT.EQ.0)THEN
         AS = AS0 / ( 1D0 + BETA0APF(NF) * AS0 * LRRAT )
*
*     Numerical solution of the evolution equation
*     with fourth-order Runge-Kutta beyond leading order
*
      ELSE
         DO K1=1,NSTEP
            XK0 = DLR * FBETA(AS              ,NF,IPT)
            XK1 = DLR * FBETA(AS + 0.5D0 * XK0,NF,IPT)
            XK2 = DLR * FBETA(AS + 0.5D0 * XK1,NF,IPT)
            XK3 = DLR * FBETA(AS +         XK2,NF,IPT)
            AS  = AS + SXTH * ( XK0 + 2D0 * XK1 + 2D0 * XK2 + XK3 )
         ENDDO
      ENDIF
*
      AS_EXACT = AS
*
      RETURN
      END
*
************************************************************************
      function as_lambda(nf,lambda2,mu2,ipt)
*
      implicit none
**
*     Input Variables
*
      integer nf,ipt
      double precision mu2,lambda2
**
*     Internal Variables
*
      double precision beta0apf,beta1apf,beta2apf
      double precision b1,b2
      double precision lo,L,lnL
**
*     Output Variables
*
      double precision as_lambda
*
      L   = dlog(mu2/lambda2)
      b1  = beta1apf(nf) / beta0apf(nf)
      b2  = beta2apf(nf) / beta0apf(nf)
      lnL = dlog(L)
      lo  = 1d0 / L / beta0apf(nf)
*
      as_lambda = lo
      if(ipt.ge.1)then
         as_lambda = as_lambda - lo**2 * b1 * lnL
      endif
      if(ipt.ge.2)then
         as_lambda = as_lambda 
     1             + lo**3 * ( b2 + b1**2 * ( lnL**2 - lnL - 1d0 ) )
      endif
*
      return
      end
*
****************************************************************************
*
*     QCD beta function.
*
****************************************************************************
      function fbeta(a,nf,ipt)
*
      implicit none
**
*     Input Variables
*
      double precision a
      integer nf,ipt
**
*     Internal Variables
*
      double precision beta0apf,beta1apf,beta2apf,beta3apf
**
*     Output Variables
*
      double precision fbeta
*
      if(ipt.eq.0)then
         fbeta = - a**2 * beta0apf(nf)
      elseif(ipt.eq.1)then
         fbeta = - a**2 * ( beta0apf(nf) + a * beta1apf(nf) )
      elseif(ipt.eq.2)then
         fbeta = - a**2 * ( beta0apf(nf)
     1           + a * ( beta1apf(nf) + a * beta2apf(nf) ) )
      elseif(ipt.ge.3)then
         fbeta = - a**2 * ( beta0apf(nf)
     1           + a * ( beta1apf(nf)
     2           + a * ( beta2apf(nf) + a * beta3apf(nf) ) ) )
      endif
*
      return
      end
*
****************************************************************************
      function beta0apf(nf)
*
      implicit none
**
*     Input Variables
*
      integer nf
**
*     Output Variables
*
      double precision beta0apf
*
      beta0apf = ( 33d0 - 2d0 * nf ) / 3d0
*
      return
      end
*
****************************************************************************
      function beta1apf(nf)
*
      implicit none
**
*     Input Variables
*
      integer nf
**
*     Output Variables
*
      double precision beta1apf
*
      beta1apf = 102d0 - 38d0 / 3d0 * nf
*
      return
      end
*
****************************************************************************
      function beta2apf(nf)
*
      implicit none
**
*     Input Variables
*
      integer nf
**
*     Output Variables
*
      double precision beta2apf
*
      beta2apf = 2857d0 / 2d0 - 5033d0 / 18d0 * nf 
     1         + 325d0 / 54d0 * nf**2
*
      return
      end
*
****************************************************************************
      function beta3apf(nf)
*
      implicit none
*
      include "../commons/consts.h"
**
*     Input Variables
*
      integer nf
**
*     Output Variables
*
      double precision beta3apf
*
      beta3apf = 149753d0 / 6d0 + 3564d0 * zeta3
     1         + ( - 1078361d0 / 162d0 - 6508d0 * zeta3 / 27d0 ) * nf
     2         + ( 50065d0 / 162d0 + 6472 * zeta3 / 81d0 ) * nf**2
     3         + 1093 / 729d0 * nf**3
*
      return
      end
*
****************************************************************************
*
*     The following routine computes the value of LambdaQCD for all the 
*     numbers of flavours and for the given perturbative order using the
*     reference value of LambdaQCD at the reference number of flavours.
*
****************************************************************************
      subroutine LambdaQCDnf
*
      implicit none
*
      include "../commons/lambda_ref_QCD.h"
      include "../commons/AlphaEvolution.h"
      include "../commons/Evs.h"
      include "../commons/Nf_FF.h"
**
*     Variables
*
      integer i
      double precision LambdaMatchUp,LambdaMatchDown,zriddr
      double precision acc,window
      external LambdaMatchUp,LambdaMatchDown
      parameter(acc=1d-12)
      parameter(window=0.3d0)
*
      if(n_ref_QCD.lt.3.or.n_ref_QCD.gt.6)then
         write(6,*) "Invalid reference number of flavours for LambdaQCD"
         write(6,*) "n_ref_QCD =",n_ref_QCD
         write(6,*) "  "
         call exit(-10)
      endif
*
      if(Evs.eq."FF")then
         if(n_ref_QCD.ne.Nf_FF)then
            write(6,*) "If the 'lambda' solution of the coupling ",
     1                 "equations is chosen"
            write(6,*) "in conjuction with the FFNS, the reference ",
     1                 "number of flavours"
            write(6,*) "for LambdaQCD must be equal to the number of ",
     1                 "flavours of the"
            write(6,*) "evolution scheme."
            write(6,*) "nf_Lambda =",n_ref_QCD
            write(6,*) "nf_FFNS =",Nf_FF
            write(6,*) "  "
            call exit(-10)
         else
            return
         endif
      endif
*
      if(n_ref_QCD.eq.3)then
         LambdaQCD(3) = lambda_ref_QCD
         do i=4,6
            LambdaQCD(i) = zriddr(LambdaMatchUp,i,LambdaQCD(i-1)*window,
     1                            LambdaQCD(i-1),acc)
         enddo
      elseif(n_ref_QCD.eq.6)then
         LambdaQCD(6) = lambda_ref_QCD
         do i=5,3,-1
            LambdaQCD(i) = zriddr(LambdaMatchDown,i,LambdaQCD(i+1),
     1                            LambdaQCD(i+1)/window,acc)
         enddo
      else
         LambdaQCD(n_ref_QCD) = lambda_ref_QCD
         do i=n_ref_QCD+1,6
            LambdaQCD(i) = zriddr(LambdaMatchUp,i,LambdaQCD(i-1)*window,
     1                            LambdaQCD(i-1),acc)
         enddo
         do i=n_ref_QCD-1,3,-1
            LambdaQCD(i) = zriddr(LambdaMatchDown,i,LambdaQCD(i+1),
     1                            LambdaQCD(i+1)/window,acc)
         enddo
      endif
*
      return
      end
*
****************************************************************************
*
*     Functions needed to find the roots to find LambdaQCD.
*
****************************************************************************
      function LambdaMatchUp(i,lambda)
*
      implicit none
*
      include "../commons/lambda_ref_QCD.h"
      include "../commons/m2th.h"
      include "../commons/kren.h"
      include "../commons/ipt.h"
      include "../commons/mass_scheme.h"
**
*     Input variables
*
      integer i
      double precision lambda
**
*     Internal variables
*
      double precision thr,au,ad,lambda2u,lambda2d
      double precision as_lambda
      double precision ln,c1,c2
**
*     Output variables
*
      double precision LambdaMatchUp
*
*     Matching condition
*
      ln = dlog( kren * k2th(i) )
*     Pole Mass
      if(mass_scheme.eq."Pole")then
         c1 = 2d0 / 3d0 * ln
         c2 = 4d0 / 9d0 * ln**2 + 38d0 / 3d0 * ln + 14d0 / 3d0
*     MSbar mass
      elseif(mass_scheme.eq."MSbar")then
         c1 = 2d0 / 3d0 * ln
         c2 = 4d0 / 9d0 * ln**2 + 22d0 / 3d0 * ln - 22d0 / 9d0
      endif
*
      thr = kren * m2th(i)
      lambda2u = lambda**2
      lambda2d = LambdaQCD(i-1)**2
      au = as_lambda(i,lambda2u,thr,ipt)
      ad = as_lambda(i-1,lambda2d,thr,ipt)
*
      if(ipt.eq.1)then
         ad = ad * ( 1d0 + c1 * ad )
      elseif(ipt.eq.2)then
         ad = ad * ( 1d0 + c1 * ad + c2 * ad**2 )
      endif
*
      LambdaMatchUp = au - ad
*
      return
      end
*
****************************************************************************
      function LambdaMatchDown(i,lambda)
*
      implicit none
*
      include "../commons/lambda_ref_QCD.h"
      include "../commons/m2th.h"
      include "../commons/kren.h"
      include "../commons/ipt.h"
      include "../commons/mass_scheme.h"
**
*     Input variables
*
      integer i
      double precision lambda
**
*     Internal variables
*
      double precision thr,au,ad,lambda2u,lambda2d
      double precision as_lambda
      double precision ln,c1,c2
**
*     Output variables
*
      double precision LambdaMatchDown
*
      ln = dlog( kren * k2th(i+1) )
*     Pole Mass
      if(mass_scheme.eq."Pole")then
         c1 = - 2d0 / 3d0 * ln
         c2 = 4d0 / 9d0 * ln**2 - 38d0 / 3d0 * ln - 14d0 / 3d0
*     MSbar mass
      elseif(mass_scheme.eq."MSbar")then
         c1 = - 2d0 / 3d0 * ln
         c2 = 4d0 / 9d0 * ln**2 - 22d0 / 3d0 * ln + 22d0 / 9d0
      endif
*
      thr = kren * m2th(i+1)
      lambda2u = LambdaQCD(i+1)**2
      lambda2d = lambda**2
      au = as_lambda(i+1,lambda2u,thr,ipt)
      ad = as_lambda(i,lambda2d,thr,ipt)
*
      if(ipt.eq.1)then
         au = au * ( 1d0 + c1 * au )
      elseif(ipt.eq.2)then
         au = au * ( 1d0 + c1 * au + c2 * au**2 )
      endif
*
      LambdaMatchDown = au - ad
*
      return
      end
*
****************************************************************************
*
*     Root finder (Numerical Recipes) slightly modified to fit the particular
*     task of finding LambdaQCD for the i-th flavour.
*
****************************************************************************
      FUNCTION zriddr(func,i,x1,x2,xacc)
      IMPLICIT NONE
      INTEGER MAXIT,i
      DOUBLE PRECISION zriddr,x1,x2,xacc,func,UNUSED
      PARAMETER (MAXIT=60,UNUSED=-1.11E30)
      EXTERNAL func
      INTEGER j
      DOUBLE PRECISION fh,fl,fm,fnew,s,xh,xl,xm,xnew
      fl=func(i,x1)
      fh=func(i,x2)
      if((fl.gt.0..and.fh.lt.0.).or.(fl.lt.0..and.fh.gt.0.))then
         xl=x1
         xh=x2
         zriddr=UNUSED
         do 11 j=1,MAXIT
            xm=0.5*(xl+xh)
            fm=func(i,xm)
            s=sqrt(fm**2-fl*fh)
            if(s.eq.0.)return
            xnew=xm+(xm-xl)*(sign(1d0,fl-fh)*fm/s)
            if (abs(xnew-zriddr).le.xacc) return
            zriddr=xnew
            fnew=func(i,zriddr)
            if (fnew.eq.0.) return
            if(sign(fm,fnew).ne.fm) then
               xl=xm
               fl=fm
               xh=zriddr
               fh=fnew
            else if(sign(fl,fnew).ne.fl) then
               xh=zriddr
               fh=fnew
            else if(sign(fh,fnew).ne.fh) then
               xl=zriddr
               fl=fnew
            else
               write(6,*) "never get here in zriddr"
               call exit(-10)
            endif
            if(abs(xh-xl).le.xacc) return
 11      continue
         write(6,*) "zriddr exceed maximum iterations"
         call exit(-10)
      else if (fl.eq.0.) then
         zriddr=x1
      else if (fh.eq.0.) then
         zriddr=x2
      else
         write(6,*) "root must be bracketed in zriddr"
         call exit(-10);
      endif
      return
      END
*
****************************************************************************
*
*     Compute and store the values of alphas right below and right below
*     the heavy quark thresholds.
*
****************************************************************************
      subroutine ThresholdAlphaQCD
*
      implicit none
*
      include "../commons/ThresholdAlphaQCD.h"
      include "../commons/kren.h"
      include "../commons/m2th.h"
**
*     Internal Variables
*
      integer inf
      double precision a_QCD
      double precision eps
      parameter(eps=1d-10)
*
      do inf=4,6
         asthUp(inf)   = a_QCD( m2th(inf) * ( 1d0 + eps ) )
         asthDown(inf) = a_QCD( m2th(inf) * ( 1d0 - eps ) )
      enddo
*
      return
      end
*
************************************************************************
      function muR2(as)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/alpha_ref_QCD.h"
      include "../commons/ThresholdAlphaQCD.h"
      include "../commons/Nf_FF.h"
      include "../commons/ipt.h"
      include "../commons/Evs.h"
      include "../commons/MaxFlavourAlpha.h"
**
*     Input Variables
*
      double precision as
**
*     Internal Variables
*
      integer nfi,nff
      integer dnf,snf
      double precision as0,muR20
      double precision integral,intBeta
      double precision eps
      parameter(eps=1d-9)
**
*     Output Variables
*
      double precision muR2
*
      as0 = ( 1d0 + eps ) * alpha_ref_QCD / 4d0 / pi
*
      muR20 = q2_ref_QCD
*
      if(Evs.eq."FF")then
         nfi = Nf_FF
         nff = Nf_FF
      elseif(Evs.eq."VF")then
         if(as.le.asthUp(6))then
            nff = 6
         elseif(as.le.asthUp(5))then
            nff = 5
         elseif(as.le.asthUp(4))then
            nff = 4
         else
            nff = 3
         endif
         if(nff.gt.nfMaxAlpha) nff = nfMaxAlpha
*
         if(as0.le.asthDown(6))then
            nfi = 6
         elseif(as0.le.asthDown(5))then
            nfi = 5
         elseif(as0.le.asthDown(4))then
            nfi = 4
         else
            nfi = 3
         endif
         if(nfi.gt.nfMaxAlpha) nfi = nfMaxAlpha
      endif
*
      integral = 0d0
 10   if(nff.eq.nfi) then
         integral = integral + intBeta(nfi,ipt,as0,as)
         muR2 = muR20 * dexp( integral )
         return
      else
         if(nff.gt.nfi)then
            dnf = 1
            snf = 1
         else
            dnf = -1
            snf = 0
         endif
*
         if(nff.gt.nfi)then
            integral = integral + intBeta(nfi,ipt,as0,asthDown(nfi+snf))
            as0 = asthUp(nfi+snf)
         else
            integral = integral + intBeta(nfi,ipt,as0,asthUp(nfi))
            as0 = asthDown(nfi)
         endif
         nfi  = nfi + dnf
         goto 10
      endif
*
      end
*
************************************************************************
*
*     Integral of the inverse of the QCD beta function.
*
************************************************************************
      function intBeta(nf,ipt,as0,as)
*
      implicit none
**
*     Input Variables
*
      integer nf,ipt
      double precision as0,as
**
*     Internal Variables
*
      double precision dgauss,eps
      double precision invBeta
      external invBeta
      parameter(eps=1d-7)

      integer nf_beta,ipt_beta
      common / BetaPar / nf_beta,ipt_beta
**
*     Output Variables
*
      double precision intBeta
*
      nf_beta  = nf
      ipt_beta = ipt
*
      intBeta = dgauss(invBeta,as0,as,eps)
*
      return
      end
*
************************************************************************
*
*     Inverse of the QCD beta function wrapped to be integrated with
*     dgauss.
*
************************************************************************
      function invBeta(as)
*
      implicit none
**
*     Input Variables
*
      double precision as
**
*     Internal Variables
*
      double precision fbeta

      integer nf_beta,ipt_beta
      common / BetaPar / nf_beta,ipt_beta
**
*     Output Variables
*
      double precision invBeta
*
      invBeta = 1d0 / fbeta(as,nf_beta,ipt_beta)
*
      return
      end
*
****************************************************************************
*
*     Set of functions to evaluate the QCD and the QED couplings simultaneuosly
*     when mixed corrections are included... it's probably not necessary.
*
****************************************************************************
      subroutine aQCDxQED_exact(ipt,nf,nl,a0,mu20,mu2,a)
*
      IMPLICIT NONE
**
*     Input Variables
*
      integer ipt,nf,nl
      double precision a0(2),mu20,mu2
**
*     Internal Variables
*
      integer nstep,istep,ia
      double precision fbetaQCDxQED
      double precision dlr
      double precision k1(2),k2(2),k3(2),k4(2)
      parameter(nstep=10)
**
*     Output Variables
*
      double precision a(2)
*
      do ia=1,2
         a(ia) = a0(ia)
      enddo
      if(mu2.eq.mu20) return
*
      dlr = dlog( mu2 / mu20 ) / nstep
*
      do istep=1,nstep
         do ia=1,2
            k1(ia) = dlr * fbetaQCDxQED(ia,ipt,nf,nl,a(1),a(2))
         enddo
         do ia=1,2
            k2(ia) = dlr * fbetaQCDxQED(ia,ipt,nf,nl,
     1           a(1) + 0.5d0 * k1(1), a(2) + 0.5d0 * k1(2))
         enddo
         do ia=1,2
            k3(ia) = dlr * fbetaQCDxQED(ia,ipt,nf,nl,
     1           a(1) + 0.5d0 * k2(1),a(2) + 0.5d0 * k2(2))
         enddo
         do ia=1,2
            k4(ia) = dlr * fbetaQCDxQED(ia,ipt,nf,nl,
     1           a(1) + k3(1),a(2) + k3(2))
         enddo
         do ia=1,2
            a(ia) = a(ia) + ( k1(ia) + 2d0 * k2(ia)
     1            + 2d0 * k3(ia) + k4(ia) ) / 6d0
         enddo
      enddo
*
      return
      end
*
****************************************************************************
*
*     QCDxQED beta function.
*
****************************************************************************
      function fbetaQCDxQED(i,ipt,nf,nl,as,a)
*
      implicit none
**
*     Input Variables
*
      integer i,ipt,nf,nl
      double precision as,a
**
*     Internal Variables
*
      double precision beta0apf,beta1apf,beta2apf
      double precision beta0qed,beta1qed
      double precision beta1qcdqed,beta1qedqcd
**
*     Output Variables
*
      double precision fbetaQCDxQED
*
      fbetaQCDxQED = 0d0
      if(ipt.eq.0)then
         if(i.eq.1)then
            fbetaQCDxQED = - as**2 * beta0apf(nf)
         elseif(i.eq.2)then
            fbetaQCDxQED = - a**2  * beta0qed(nf,nl)
         endif
      elseif(ipt.eq.1)then
         if(i.eq.1)then
            fbetaQCDxQED = - as**2 * ( beta0apf(nf)
     1           + as * beta1apf(nf) + a * beta1qcdqed(nf) )
         elseif(i.eq.2)then
            fbetaQCDxQED = - a**2 * ( beta0qed(nf,nl)
     1        + a * beta1qed(nf,nl) + as * beta1qedqcd(nf) )
         endif
      elseif(ipt.eq.2)then
         if(i.eq.1)then
            fbetaQCDxQED = - as**2 * ( beta0apf(nf)
     1           + as * beta1apf(nf) + a * beta1qcdqed(nf)
     2           + as**2 * beta2apf(nf) )
         elseif(i.eq.2)then
            fbetaQCDxQED = - a**2 * ( beta0qed(nf,nl)
     1        + a * beta1qed(nf,nl) + as * beta1qedqcd(nf) )
         endif
      endif
*
      return
      end
