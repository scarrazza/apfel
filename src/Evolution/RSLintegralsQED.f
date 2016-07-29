************************************************************************
*
*     RSLintegralsQED.f:
*
*     This routine evaluates and dumps on the common integrals.h once and
*     for all for a given number of active flavours nf and for the pair 
*     (alpha,beta) the integral of the sum of Regular and Singular part 
*     plus the local part of all the splitting functions for all the needed
*     orders in QED (including alphas corrections).
*
*     The index kk runs like that:
*
*     kk  =  1   2   3   4   5   6   7   8   9   10  11  12  13  14
*            nsp nsm gg  gq  gD  qg  qq  qD  Dg  Dq  DD  LL  gL  Lg
* 
************************************************************************
      subroutine RSLintegralsQED(nf,nl,beta,alpha)
*
      implicit none
*
      include "../commons/ColorFactors.h"
      include "../commons/wrap.h"
      include "../commons/grid.h"
      include "../commons/integrals.h"
      include "../commons/LeptEvol.h"
      include "../commons/ipt.h"
**
*     Input Variables
*
      integer nf,nl,beta,alpha
**
*     Internal Variables
*
      integer bound
      integer nfup,nfdn,nc,nk
      integer jpt
      double precision fL
      double precision PL(0:2),integ(0:2),cp(0:2)
      double precision X0NSC
      double precision X1NSC_ASA,X1GAMGAMC_ASA
      double precision X1GAMGAMC_AA,REM_X1NSC_AA
      double precision dgauss,a,b,eps(0:1)
      double precision integrandsQED
      double precision e2u,e2d,e2sig,de2,etap,etam
      double precision e4u,e4d,e4sig,de4
      double precision ns0L,ns0RL,qg0R,gq0R
      double precision ns1L,gmgm1L
      double precision ns2L,gmgm2L
      external integrandsQED
      data eps / 1d-6, 1d-5 /
      parameter(nc = 3)
      parameter(e2u =  4d0 / 9d0)
      parameter(e2d =  1d0 / 9d0)
      parameter(e4u = 16d0 / 81d0)
      parameter(e4d =  1d0 / 81d0)
*
      jpt = 0
      if(ipt.ge.1) jpt = 2
*
*     Initialize Integrals
*
      do k=1,14
         do wipt=0,jpt
            SQ(igrid,nf,nl,k,wipt,beta,alpha) = 0d0
         enddo
      enddo
*
*     Adjustment od the bounds of the integrals
*
      if(alpha.lt.beta)then
         return
      else
         bound = alpha-inter_degree(igrid)
         if(alpha.lt.inter_degree(igrid)) bound = 0
         a = max(xg(igrid,beta),xg(igrid,beta)/xg(igrid,alpha+1))
         b = min(1d0,xg(igrid,beta)/xg(igrid,bound))
      endif
*
      fL = 0d0
      if(alpha.eq.beta) fL = 1d0
*
*     Couplings
*
      if(nf.eq.3)then
         nfup = 1
         nfdn = 2
      elseif(nf.eq.4)then
         nfup = 2
         nfdn = 2
      elseif(nf.eq.5)then
         nfup = 2
         nfdn = 3
      elseif(nf.eq.6)then
         nfup = 3
         nfdn = 3
      endif
*
      e2sig = nc * ( nfup * e2u + nfdn * e2d )
      de2   = nc * ( nfup * e2u - nfdn * e2d )
      if(LeptEvol) e2sig = e2sig + nl
*
      etap = ( e2u + e2d ) / 2d0
      etam = ( e2u - e2d ) / 2d0
*
*     Variables needed for wrapping the integrand functions
*
      walpha = alpha
      wbeta  = beta
*
*     Precompute integrals
*
      wipt = 0
      k = 1
      ns0RL = dgauss(integrandsQED,a,b,eps(0))
      k = 6
      qg0R  = dgauss(integrandsQED,a,b,eps(0))
      k = 4
      gq0R  = dgauss(integrandsQED,a,b,eps(0))
      ns0L  = X0NSC(a)
*     O(alpha_s alpha) contributions
      if(jpt.ge.1)then
         e4sig = nc * ( nfup * e4u + nfdn * e4d )
         de4   = nc * ( nfup * e4u - nfdn * e4d )
         if(LeptEvol) e4sig = e4sig + nl
         wipt = 1

         ns1L   = X1NSC_ASA()
         gmgm1L = X1GAMGAMC_ASA()
      endif
*     O(alpha^2) contributions
      if(jpt.ge.1)then
         wipt = 2
         
         ns2L   = REM_X1NSC_AA(a)
         gmgm2L = X1GAMGAMC_AA()
      endif
*
      nk = 11
      if(LeptEvol) nk = 14
*
      do k=1,nk
*
*     LO
*
*     Non-singlet Plus
         if(k.eq.1)then
            cp(0) = e2u
            PL(0) = ns0L / CF
            integ(0) = ns0RL
*     Non-singlet Minus
         elseif(k.eq.2)then
            cp(0) = e2d
            PL(0) = ns0L / CF
            integ(0) = ns0RL
*    Gamma-Gamma
         elseif(k.eq.3)then
            cp(0) = e2sig
            PL(0) = - 4d0 / 3d0
            integ(0) = 0d0
*    Gamma-Quark
         elseif(k.eq.4)then
            cp(0) = etap
            PL(0) = 0d0
            integ(0) = gq0R
*    Gamma-Delta
         elseif(k.eq.5)then
            cp(0) = etam
            PL(0) = 0d0
            integ(0) = gq0R
*    Quark-Gamma
         elseif(k.eq.6)then
            cp(0) = 2d0 * e2sig
            PL(0) = 0d0
            integ(0) = qg0R
*    Quark-Quark
         elseif(k.eq.7)then
            cp(0) = etap
            PL(0) = ns0L / CF
            integ(0) = ns0RL
*    Quark-Delta
         elseif(k.eq.8)then
            cp(0) = etam
            PL(0) = ns0L / CF
            integ(0) = ns0RL
*    Delta-Gamma
         elseif(k.eq.9)then
            cp(0) = 2d0 * de2
            PL(0) = 0d0
            integ(0) = qg0R
*    Delta-Quark
         elseif(k.eq.10)then
            cp(0) = etam
            PL(0) = ns0L / CF
            integ(0) = ns0RL
*    Delta-Delta
         elseif(k.eq.11)then
            cp(0) = etap
            PL(0) = ns0L / CF
            integ(0) = ns0RL
*     Lepton-Lepton
         elseif(k.eq.12)then
            cp(0) = 1d0
            PL(0) = ns0L / CF
            integ(0) = ns0RL
*    Gamma-Lepton
         elseif(k.eq.13)then
            cp(0) = 1d0
            PL(0) = 0d0
            integ(0) = gq0R
*    Lepton-Gamma
         elseif(k.eq.14)then
            cp(0) = 1d0
            PL(0) = 0d0
            integ(0) = 2d0 * nl * qg0R
         endif
*
*     NLO (i.e. O(alpha_s alpha))
*
         if(jpt.ge.1)then


         endif
*
*     NNLO (i.e. O(alpha^2))
*
         if(jpt.ge.2)then

         endif
*
*     Integrals
*
         do wipt=0,jpt
            SQ(igrid,nf,nl,k,wipt,beta,alpha) = cp(wipt)
     1           * ( integ(wipt) + PL(wipt) * fL )
         enddo
      enddo
*
      return
      end
