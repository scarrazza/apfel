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
*     kk  =  1    2    3    4
*            nspu nsmu nspd nsmd
*
*            5    6    7    8
*            gg   ggm  gS   gD
*
*            9    10   11   12
*            gmg  gmgm gmS  gmD
*
*            13   14   15   16
*            Sg   Sgm  SS   SD
*
*            17   18   19   20
*            Dg   Dgm  DS   DD
*
*            21   22
*            VV   VDV (DVDV = VV, DVV = VDV)
*
*            23   24   25
*            LL   gmL  Lgm
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
      double precision dgauss,a,b,eps(0:2)
      double precision integrandsQED
      double precision e2u,e2d,e2sig,de2,etap,etam
      double precision e4u,e4d,e4sig,de4
      double precision ns0L,ns0RS,qg0R,gq0R
      double precision ns1L,gmgm1L,nsp1RS,nsm1RS,qg1R,gq1R,ggm1R
      double precision ns2L,gmgm2L,ns2RS,gq2R,ps2R
      external integrandsQED
      data eps / 1d-6, 1d-5, 1d-5 /
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
      ns0RS = dgauss(integrandsQED,a,b,eps(0))
      k = 14
      qg0R  = dgauss(integrandsQED,a,b,eps(0))
      k = 11
      gq0R  = dgauss(integrandsQED,a,b,eps(0))
      ns0L  = X0NSC(a)
*     O(alpha_s alpha) contributions
      if(jpt.ge.1)then
         wipt = 1
*     Local pieces
         ns1L   = X1NSC_ASA()
         gmgm1L = X1GAMGAMC_ASA()
*     NS+ up, NS+ down, Quark-Quark, Quark-Delta, Delta-Quark, Delta-Delta
         k = 1
         nsp1RS = dgauss(integrandsQED,a,b,eps(1))
*     NS- up, NS- down, Valence-Valence, Valence-D_V
         k = 3
         nsm1RS = dgauss(integrandsQED,a,b,eps(1))
*     Quark-Gluon, Delta-Gluon, Quark-Gamma, Delta-Gamma
         k = 13
         qg1R = dgauss(integrandsQED,a,b,eps(1))
*     Gluon-Quark, Gluon-Delta, Gamma-Quark, Gamma-Delta
         k = 7
         gq1R = dgauss(integrandsQED,a,b,eps(1))
*     Gluon-Gamma, Gamma-Gluon
         k = 6
         ggm1R = dgauss(integrandsQED,a,b,eps(1))
      endif
*     O(alpha^2) contributions
      if(jpt.ge.2)then
         e4sig = nc * ( nfup * e4u + nfdn * e4d )
         de4   = nc * ( nfup * e4u - nfdn * e4d )
         if(LeptEvol) e4sig = e4sig + nl
*
         wipt = 2
*     Local pieces
         ns2L   = REM_X1NSC_AA(a)
         gmgm2L = X1GAMGAMC_AA()
*     Remainder NS
         k = 1
         ns2RS  = dgauss(integrandsQED,a,b,eps(2))
*     Remainder Gamma-Quark, Gamma-Delta
         k = 11
         gq2R = dgauss(integrandsQED,a,b,eps(2))
*     Pure singlet
         k = 15
         ps2R = dgauss(integrandsQED,a,b,eps(2))
      endif
*
      nk = 22
      if(LeptEvol) nk = 25
*
      do k=1,nk
*
*     LO
*
         cp(0) = 0d0
         PL(0) = 0d0
         integ(0) = 0d0
*     NS+ up
         if(k.eq.1)then
            cp(0) = e2u
            PL(0) = ns0L / CF
            integ(0) = ns0RS
*     NS+ down
         elseif(k.eq.2)then
            cp(0) = e2d
            PL(0) = ns0L / CF
            integ(0) = ns0RS
*     NS- up
         elseif(k.eq.3)then
            cp(0) = e2u
            PL(0) = ns0L / CF
            integ(0) = ns0RS
*     NS- down
         elseif(k.eq.4)then
            cp(0) = e2d
            PL(0) = ns0L / CF
            integ(0) = ns0RS
*    Gamma-Gamma
         elseif(k.eq.10)then
            cp(0) = e2sig
            PL(0) = - 4d0 / 3d0
            integ(0) = 0d0
*    Gamma-Quark
         elseif(k.eq.11)then
            cp(0) = etap
            PL(0) = 0d0
            integ(0) = gq0R
*    Gamma-Delta
         elseif(k.eq.12)then
            cp(0) = etam
            PL(0) = 0d0
            integ(0) = gq0R
*    Quark-Gamma
         elseif(k.eq.14)then
            cp(0) = 2d0 * e2sig
            PL(0) = 0d0
            integ(0) = qg0R
*    Quark-Quark
         elseif(k.eq.15)then
            cp(0) = etap
            PL(0) = ns0L / CF
            integ(0) = ns0RS
*    Quark-Delta
         elseif(k.eq.16)then
            cp(0) = etam
            PL(0) = ns0L / CF
            integ(0) = ns0RS
*    Delta-Gamma
         elseif(k.eq.18)then
            cp(0) = 2d0 * de2
            PL(0) = 0d0
            integ(0) = qg0R
*    Delta-Quark
         elseif(k.eq.19)then
            cp(0) = etam
            PL(0) = ns0L / CF
            integ(0) = ns0RS
*    Delta-Delta
         elseif(k.eq.20)then
            cp(0) = etap
            PL(0) = ns0L / CF
            integ(0) = ns0RS
*    Valence-Valence (=D_V-D_V)
         elseif(k.eq.21)then
            cp(0) = etap
            PL(0) = ns0L / CF
            integ(0) = ns0RS
*    Valence-D_V (=D_V-Valence)
         elseif(k.eq.22)then
            cp(0) = etam
            PL(0) = ns0L / CF
            integ(0) = ns0RS
*     Lepton-Lepton
         elseif(k.eq.23)then
            cp(0) = 1d0
            PL(0) = ns0L / CF
            integ(0) = ns0RS
*    Gamma-Lepton
         elseif(k.eq.24)then
            cp(0) = 1d0
            PL(0) = 0d0
            integ(0) = gq0R
*    Lepton-Gamma
         elseif(k.eq.25)then
            cp(0) = 1d0
            PL(0) = 0d0
            integ(0) = 2d0 * nl * qg0R
         endif
*
*     NLO (i.e. O(alpha_s alpha))
*
         if(jpt.ge.1)then
            cp(1) = 0d0
            PL(1) = 0d0
            integ(1) = 0d0
*     NS+ up
            if(k.eq.1)then
               cp(1) = e2u
               PL(1) = ns1L
               integ(1) = nsp1RS
*     NS+ down
            elseif(k.eq.2)then
               cp(1) = e2d
               PL(1) = ns1L
               integ(1) = nsp1RS
*     NS- up
            elseif(k.eq.3)then
               cp(1) = e2u
               PL(1) = ns1L
               integ(1) = nsm1RS
*     NS- down
            elseif(k.eq.4)then
               cp(1) = e2d
               PL(1) = ns1L
               integ(1) = nsm1RS
*     Gluon-Gluon
            elseif(k.eq.5)then
               cp(1) = e2sig
               PL(1) = gmgm1L * TR / CF
               integ(1) = 0d0
*     Gluon-Gamma
            elseif(k.eq.6)then
               cp(1) = e2sig
               PL(1) = 0d0
               integ(1) = ggm1R
*     Gluon-Quark
            elseif(k.eq.7)then
               cp(1) = etap
               PL(1) = 0d0
               integ(1) = gq1R
*     Gluon-Delta
            elseif(k.eq.8)then
               cp(1) = etam
               PL(1) = 0d0
               integ(1) = gq1R
*     Gamma-Gluon
            elseif(k.eq.9)then
               cp(1) = e2sig
               PL(1) = 0d0
               integ(1) = ggm1R * TR / CF
*     Gamma-Gamma
            elseif(k.eq.10)then
               cp(1) = e2sig
               PL(1) = gmgm1L
               integ(1) = 0d0
*     Gamma-Quark
            elseif(k.eq.11)then
               cp(1) = etap
               PL(1) = 0d0
               integ(1) = gq1R
*     Gamma-Delta
            elseif(k.eq.12)then
               cp(1) = etam
               PL(1) = 0d0
               integ(1) = gq1R
*     Quark-Gluon
            elseif(k.eq.13)then
               cp(1) = 2d0 * e2sig
               PL(1) = 0d0
               integ(1) = qg1R * TR / CF
*     Quark-Gamma
            elseif(k.eq.14)then
               cp(1) = 2d0 * e2sig
               PL(1) = 0d0
               integ(1) = qg1R
*     Quark-Quark
            elseif(k.eq.15)then
               cp(1) = etap
               PL(1) = ns1L
               integ(1) = nsp1RS
*     Quark-Delta
            elseif(k.eq.16)then
               cp(1) = etam
               PL(1) = ns1L
               integ(1) = nsp1RS
*     Delta-Gluon
            elseif(k.eq.17)then
               cp(1) = 2d0 * de2
               PL(1) = 0d0
               integ(1) = qg1R * TR / CF
*     Delta-Gamma
            elseif(k.eq.18)then
               cp(1) = 2d0 * de2
               PL(1) = 0d0
               integ(1) = qg1R
*     Delta-Quark
            elseif(k.eq.19)then
               cp(1) = etam
               PL(1) = ns1L
               integ(1) = nsp1RS
*     Delta-Delta
            elseif(k.eq.20)then
               cp(1) = etap
               PL(1) = ns1L
               integ(1) = nsp1RS
*     Valence-Valence (=D_V-D_V)
            elseif(k.eq.21)then
               cp(1) = etap
               PL(1) = ns1L
               integ(1) = nsm1RS
*     Valence-D_V (=D_V-Valence)
            elseif(k.eq.22)then
               cp(1) = etam
               PL(1) = ns1L
               integ(1) = nsm1RS
            endif
         endif
*
*     NNLO (i.e. O(alpha^2))
*
         if(jpt.ge.2)then
            cp(2) = 0d0
            PL(2) = 0d0
            integ(2) = 0d0
*     NS+ up
            if(k.eq.1)then
               cp(2) = e4u
               PL(2) = ns1L / 2d0 / CF + e2sig / e2u * ns2L
               integ(2) = nsp1RS / 2d0 / CF + e4sig / e2u * ns2RS
*     NS+ down
            elseif(k.eq.2)then
               cp(2) = e4d
               PL(2) = ns1L / 2d0 / CF + e2sig / e2d * ns2L
               integ(2) = nsp1RS / 2d0 / CF + e4sig / e2d * ns2RS
*     NS- up
            elseif(k.eq.3)then
               cp(2) = e4u
               PL(2) = ns1L / 2d0 / CF + e2sig / e2u * ns2L
               integ(2) = nsm1RS / 2d0 / CF + e4sig / e2u * ns2RS
*     NS- down
            elseif(k.eq.4)then
               cp(2) = e4d
               PL(2) = ns1L / 2d0 / CF + e2sig / e2d * ns2L
               integ(2) = nsm1RS / 2d0 / CF + e4sig / e2d * ns2RS
*     Gamma-Gamma
            elseif(k.eq.10)then
               cp(2) = e4sig
               PL(2) = gmgm1L
               integ(2) = ggm1R / CF
*     Gamma-Quark
            elseif(k.eq.11)then
               cp(2) = 1d0
               PL(2) = 0d0
               integ(2) = gq1R / CF * ( e4u + e4d ) / 2d0
     1                  + e2sig * gq2R * ( e2u + e2d ) / 2d0
*     Gamma-Delta
            elseif(k.eq.12)then
               cp(2) = 1d0
               PL(2) = 0d0
               integ(2) = gq1R / CF * ( e4u - e4d ) / 2d0
     1                  + e2sig * gq2R * ( e2u - e2d ) / 2d0
*     Quark-Gamma
            elseif(k.eq.14)then
               cp(2) = 2d0 * e4sig
               PL(2) = 0d0
               integ(2) = qg1R / CF
*     Quark-Quark
            elseif(k.eq.15)then
               cp(2) = 1d0
               PL(2) = ns1L / 2d0 / CF * ( e4u + e4d ) / 2d0
     1               + e2sig * ns2L * ( e2u + e2d ) / 2d0
               integ(2) = nsp1RS / 2d0 / CF * ( e4u + e4d ) / 2d0
     1                  + e2sig * ns2RS * ( e2u + e2d ) / 2d0
     2                  + etap * e2sig * ps2R
*     Quark-Delta
            elseif(k.eq.16)then
               cp(2) = 1d0
               PL(2) = ns1L / 2d0 / CF * ( e4u - e4d ) / 2d0
     1               + e2sig * ns2L * ( e2u - e2d ) / 2d0
               integ(2) = nsp1RS / 2d0 / CF * ( e4u - e4d ) / 2d0
     1                  + e2sig * ns2RS * ( e2u - e2d ) / 2d0
     2                  + etam * e2sig * ps2R
*     Delta-Gamma
            elseif(k.eq.18)then
               cp(2) = 2d0 * de4
               PL(2) = 0d0
               integ(2) = qg1R / CF
*     Delta-Quark
            elseif(k.eq.19)then
               cp(2) = 1d0
               PL(2) = ns1L / 2d0 / CF * ( e4u - e4d ) / 2d0
     1               + e2sig * ns2L * ( e2u - e2d ) / 2d0
               integ(2) = nsp1RS / 2d0 / CF * ( e4u - e4d ) / 2d0
     1                  + e2sig * ns2RS * ( e2u - e2d ) / 2d0
     2                  + etam * de2 * ps2R
*     Delta-Delta
            elseif(k.eq.20)then
               cp(2) = 1d0
               PL(2) = ns1L / 2d0 / CF * ( e4u + e4d ) / 2d0
     1               + e2sig * ns2L * ( e2u + e2d ) / 2d0
               integ(2) = nsp1RS / 2d0 / CF * ( e4u + e4d ) / 2d0
     1                  + e2sig * ns2RS * ( e2u + e2d ) / 2d0
     2                  + etap * de2 * ps2R
*     Valence-Valence (=D_V-D_V)
            elseif(k.eq.21)then
               cp(2) = 1d0
               PL(2) = ns1L / 2d0 / CF * ( e4u + e4d ) / 2d0
     1               + e2sig * ns2L * ( e2u + e2d ) / 2d0
               integ(2) = nsm1RS / 2d0 / CF * ( e4u + e4d ) / 2d0
     1                  + e2sig * ns2RS * ( e2u + e2d ) / 2d0
*     Valence-D_V (=D_V-Valence)
            elseif(k.eq.22)then
               cp(2) = 1d0
               PL(2) = ns1L / 2d0 / CF * ( e4u - e4d ) / 2d0
     1               + e2sig * ns2L * ( e2u - e2d ) / 2d0
               integ(2) = nsm1RS / 2d0 / CF * ( e4u - e4d ) / 2d0
     1                  + e2sig * ns2RS * ( e2u - e2d ) / 2d0
            endif
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
