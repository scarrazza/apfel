************************************************************************
*
*     RSLintegralsQED.f:
*
*     This routine evaluates and dump on the common integrals.h once and 
*     for all for a given number of active flavours nf and for the pair 
*     (alpha,beta) the integral of the sum of Regular and Singular part 
*     plus the local part of all the splitting functions for all the needed
*     orders in QED.
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
      include "../commons/wrap.h"
      include "../commons/grid.h"
      include "../commons/integrals.h"
**
*     Input Variables
*
      integer nf,nl,beta,alpha
**
*     Internal Variables
*
      integer bound
      integer nfup,nfdn,nc
      double precision PL(14),fL
      double precision X0NSC
      double precision dgauss,a,b,eps
      double precision integrandsQED
      double precision e2u,e2d,e2sig,fnf,etap,etam,thetap,thetam,CF
      double precision nsL,nsRL,qgR,gqR,integ,cp(14)
      external integrandsQED
      parameter(nc=3)
      parameter(e2u=4d0/9d0)
      parameter(e2d=1d0/9d0)
      parameter(CF=4d0/3d0)
      parameter(eps=1d-6)
*
*     Initialize Integrals
*
      do k=1,14
         SQ(igrid,nf,nl,k,beta,alpha) = 0d0
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
*
      etap = ( e2u + e2d ) / 2d0
      etam = ( e2u - e2d ) / 2d0
*
      fnf = dble( nfup - nfdn ) / dble(nf)
      thetap = 2d0 * nc * nf * ( fnf * etap + etam )
      thetam = 2d0 * nc * nf * ( etap + fnf * etam )
*
*     Variables needed for wrapping the integrand functions
*
      walpha = alpha
      wbeta  = beta
*
*     Precompute integrals
*
      k = 1
      nsRL = dgauss(integrandsQED,a,b,eps)
      k = 6
      qgR  = dgauss(integrandsQED,a,b,eps)
      k = 4
      gqR  = dgauss(integrandsQED,a,b,eps)
      nsL  = X0NSC(a)
*
      do k=1,14
*
*     Local contributions
*
*     Non-singlet Plus
         if(k.eq.1)then
            cp(k) = e2u
            PL(k) = nsL / CF
            integ = nsRL
*     Non-singlet Minus
         elseif(k.eq.2)then
            cp(k) = e2d
            PL(k) = nsL / CF
            integ = nsRL
*    Gamma-Gamma
         elseif(k.eq.3)then
            cp(k) = e2sig
            PL(k) = - 4d0 / 3d0
            integ = 0d0
*    Gamma-Quark
         elseif(k.eq.4)then
            cp(k) = etap
            PL(k) = 0d0
            integ = gqR
*    Gamma-Delta
         elseif(k.eq.5)then
            cp(k) = etam
            PL(k) = 0d0
            integ = gqR
*    Quark-Gamma
         elseif(k.eq.6)then
            cp(k) = thetam
            PL(k) = 0d0
            integ = qgR
*    Quark-Quark
         elseif(k.eq.7)then
            cp(k) = etap
            PL(k) = nsL / CF
            integ = nsRL
*    Quark-Delta
         elseif(k.eq.8)then
            cp(k) = etam
            PL(k) = nsL / CF
            integ = nsRL
*    Delta-Gamma
         elseif(k.eq.9)then
            cp(k) = thetap
            PL(k) = 0d0
            integ = qgR
*    Delta-Quark
         elseif(k.eq.10)then
            cp(k) = etam
            PL(k) = nsL / CF
            integ = nsRL
*    Delta-Delta
         elseif(k.eq.11)then
            cp(k) = etap
            PL(k) = nsL / CF
            integ = nsRL
*     Lepton-Lepton
         elseif(k.eq.12)then
            cp(k) = 1d0
            PL(k) = nsL / CF
            integ = nsRL
*    Gamma-Lepton
         elseif(k.eq.13)then
            cp(k) = 1d0
            PL(k) = 0d0
            integ = gqR
*    Lepton-Gamma
         elseif(k.eq.14)then
            cp(k) = 1d0
            PL(k) = 0d0
            integ = 2d0 * nl * qgR
         endif
*
*     Integrals
*
         SQ(igrid,nf,nl,k,beta,alpha) = cp(k) * ( integ + PL(k) * fL )
      enddo
*
      return
      end
