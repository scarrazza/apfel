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
      asr0 = alpha_ref_QCD / 4d0 / pi
*
      mur20 = q2_ref_QCD
      mur2  = kren * mu2F
      do i=4,6
         mur2th(i) = kren * m2th(i)
c         mur2th(i) = m2th(i)
      enddo
      do i=3,6
         lambda2(i) = LambdaQCD(i)**2d0
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
*     c1 and c2 are the same coefficients used in eq. (2.42) of hep-ph/0408244 and 
*     obtained in eq. (10) of hep-ph/9706430. In the following they are divided by 
*     (4*pi) and (4*pi)^2 respectively to match the notations. Note that in terms 
*     of the MSbar mass this coefficients change.
*
      kappa = kren         ! mu_R / mu_F
c      kappa = 1d0          ! mu_R / mu_F
      ln = dlog(kappa)
*     Pole Mass
      if(mass_scheme.eq."Pole")then
         if(nff.gt.nfi)then
            c1 = 2d0 / 3d0 * ln
            c2 = 4d0 / 9d0 * ln**2d0 + 38d0 / 3d0 * ln + 14d0 / 3d0
         elseif(nff.lt.nfi)then
            c1 = - 2d0 / 3d0 * ln
            c2 = 4d0 / 9d0 * ln**2d0 - 38d0 / 3d0 * ln - 14d0 / 3d0
         endif
*     MSbar mass
      elseif(mass_scheme.eq."MSbar")then
         if(nff.gt.nfi)then
            c1 = 0d0
            c2 = - 22d0 / 9d0
         elseif(nff.lt.nfi)then
            c1 = 0d0
            c2 = 22d0 / 9d0
         endif
      endif
*
 10   if(nff.eq.nfi) then
         a_QCD = as_exact(nfi,mur20,asr0,mur2,ipt)
c         a_QCD = as_expanded(nfi,mur20,asr0,mur2,ipt)
c         a_QCD = as_lambda(nfi,lambda2(nfi),mur2,ipt)
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
         asi = as_exact(nfi,mur20,asr0,mur2th(nfi+snf),ipt)
c         asi = as_expanded(nfi,mur20,asr0,mur2th(nfi+snf),ipt)
c         asi = as_lambda(nfi,lambda2(nfi),mur2th(nfi+snf),ipt)
*
*     NLO and NNLO threshold matchings
*
         if(ipt.eq.1)then
            asi = asi * ( 1d0 + c1 * asi )
         elseif(ipt.eq.2)then
            asi = asi * ( 1d0 + c1 * asi + c2 * asi**2d0 )
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
      double precision beta0,beta1,beta2,b1,b2
**
*     Output Variables
*
      double precision as_expanded
*
      b1    = beta1(nf) / beta0(nf)
      b2    = beta2(nf) / beta0(nf)
*
      asi = as0
      t = log(mu2/mu20)
      den = 1d0 + beta0(nf) * asi * t
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
     1      + ( alo * ( alo - asi ) * ( b2 - b1**2d0 )
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
*
      include "../commons/consts.h"
**
*     Input Variables
*
      INTEGER NF,IPT
      DOUBLE PRECISION MU2
      DOUBLE PRECISION AS0,MU20
**
*     Inernal Variables
*
      INTEGER NFMIN,NFMAX,NSTEP,K1
      DOUBLE PRECISION AS
      DOUBLE PRECISION BETA0,FBETA
      DOUBLE PRECISION DLR,LRRAT,SXTH
      DOUBLE PRECISION XK0,XK1,XK2,XK3

      PARAMETER(NFMIN=3,NFMAX=6)
      PARAMETER(NSTEP=50)
      PARAMETER(SXTH=0.166666666666666D0 )
**
*     Output Variables
*
      DOUBLE PRECISION AS_EXACT
*
*     ..The beta functions FBETAn at N^nLO for n = 1, 2, and 3
*
      AS = AS0
      LRRAT = DLOG (MU2/MU20)
      DLR = LRRAT / NSTEP
*     
*     ..Solution of the evolution equation depending on  NAORD
*   (fourth-order Runge-Kutta beyond the leading order)
*     
      IF(IPT.EQ.0)THEN
         AS = AS0 / ( 1D0 + BETA0(NF) * AS0 * LRRAT )
      ELSEIF(IPT.EQ.1)THEN
         DO 2 K1=1,NSTEP
            XK0 = DLR * FBETA(AS,NF,IPT)
            XK1 = DLR * FBETA(AS + 0.5d0 * XK0,NF,IPT)
            XK2 = DLR * FBETA(AS + 0.5d0 * XK1,NF,IPT)
            XK3 = DLR * FBETA(AS + XK2,NF,IPT)
            AS = AS + SXTH * ( XK0 + 2D0 * XK1 + 2d0 * XK2 + XK3 )
 2        CONTINUE
      ELSEIF(IPT.EQ.2)THEN
         DO 3 K1 = 1, NSTEP
            XK0 = DLR * FBETA(AS,NF,IPT)
            XK1 = DLR * FBETA(AS + 0.5d0 * XK0,NF,IPT)
            XK2 = DLR * FBETA(AS + 0.5d0 * XK1,NF,IPT)
            XK3 = DLR * FBETA(AS + XK2,NF,IPT)
            AS = AS + SXTH * ( XK0 + 2D0 * XK1 + 2D0 * XK2 + XK3 )
 3       CONTINUE
      ENDIF
*
      AS_EXACT = AS
*
      RETURN
      END
*
************************************************************************
*
*     The following function computes alphas / 4*pi given the number of
*     active favours nf, the value of LambdaQCD^2, the final scale and
*     the perturbative order.
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
      double precision beta0,beta1,beta2
      double precision b1,b2
      double precision lo,L,lnL
**
*     Output Variables
*
      double precision as_lambda
*
      L   = dlog(mu2/lambda2)
      b1  = beta1(nf) / beta0(nf)
      b2  = beta2(nf) / beta0(nf)
      lnL = dlog(L)
      lo  = 1d0 / L / beta0(nf)
*
      as_lambda = lo
      if(ipt.ge.1)then
         as_lambda = as_lambda - lo**2d0 * b1 * lnL
      endif
      if(ipt.ge.2)then
         as_lambda = as_lambda 
     1             + lo**3d0 * ( b2 + b1**2 * ( lnL**2d0 - lnL - 1d0 ) )
      endif
*
      return
      end
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
      double precision beta0,beta1,beta2
**
*     Output Variables
*
      double precision fbeta
*
      if(ipt.eq.1)then
         fbeta = - a**2d0 * ( beta0(nf) + a * beta1(nf) )
      elseif(ipt.eq.2)then
         fbeta = - a**2d0 * ( beta0(nf) 
     1           + a * ( beta1(nf) + a * beta2(nf) ) )
      endif
*
      return
      end
*
****************************************************************************
      function beta0(nf)
*
      implicit none
**
*     Input Variables
*
      integer nf
**
*     Output Variables
*
      double precision beta0
*
      beta0 = ( 33d0 - 2d0 * nf ) / 3d0
*
      return
      end
*
****************************************************************************
      function beta1(nf)
*
      implicit none
**
*     Input Variables
*
      integer nf
**
*     Output Variables
*
      double precision beta1
*
      beta1 = 102d0 - 38d0 / 3d0 * nf
*
      return
      end
*
****************************************************************************
      function beta2(nf)
*
      implicit none
**
*     Input Variables
*
      integer nf
**
*     Output Variables
*
      double precision beta2
*
      beta2 = 2857d0 / 2d0 - 5033d0 / 18d0 * nf + 325d0 / 54d0 * nf**2d0
*
      return
      end
