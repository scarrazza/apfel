************************************************************************
*
*     RSLintegralsSIA.f:
*
*     This routine evaluates and dump on the common integralsDIS.h once and 
*     for all for the pair (alpha,beta) the integral of the sum of Regular
*     and Singular part plus the Local part of all the coefficient functions
*     for all the needed orders in SIA.
*     The index sf runs like that:
*
*        sf  Coefficient Function
*     --------------------------
*        1           C2
*        2           CL
*        3           C3

*
*     while the index kk runs like that:
*
*        kk     combination
*     --------------------------
*        1         gluon  
*        2      pure-singlet
*        3    non-singlet-plus
*        4    non-singlet-minus
*
************************************************************************
      subroutine RSLintegralsSIA(beta,alpha)
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/wrapDIS.h"
      include "../commons/grid.h"
      include "../commons/coeffhqmellin.h"
      include "../commons/integralsDIS.h"
      include "../commons/Nf_FF.h"
      include "../commons/ColorFactors.h"
**
*     Input Variables
*
      integer beta,alpha
**
*     Internal Variables
*
      integer bound,inf
      double precision C2L(4,0:2),CLL(4,0:2),C3L(4,0:2)
      double precision fL
      double precision dgauss,a,b,eps(2)
      double precision integrandsSIAzm

      double precision integC2(0:2),integCL(0:2),integC3(0:2)
      double precision C2ns1RS,C2ns1L,C2g1R,CLns1RS,CLg1R,C3ns1RS,C3ns1L
      double precision CLns1L,C3g1R
c      double precision C2g2R,C2g2L,CLg2R,C2ps2R,CLps2R,C3ps2R
c      double precision C2ns2RS,CLns2RS,C3ns2RS
c      double precision C2ns2L,C3ns2L
c      double precision C2nsp2RS,CLnsp2RS,C3nsp2RS
c      double precision C2nsp2L,CLnsp2L,C3nsp2L
c      double precision C2nsm2RS,CLnsm2RS,C3nsm2RS
c      double precision C2nsm2L,CLnsm2L,C3nsm2L
      double precision C2NS1TC,C3NS1TC
c      double precision C2NSP2TC,CLNSP2TC,C3NSP2TC,C2NSM2TC,CLNSM2TC,C3NSM2TC
      external integrandsSIAzm
c      data eps / 5d-8, 1d-3 /
c      data eps / 5d-8, 1d-5 /
      data eps / 1d-6, 1d-3 /
c      data eps / 1d-7, 1d-5 /
*
*     Adjustment of the bounds of the integrals
*
      if(alpha.lt.beta) return
*
      bound = alpha-inter_degree(igrid)
      if(alpha.lt.inter_degree(igrid)) bound = 0
      a = max(xg(igrid,beta),xg(igrid,beta)/xg(igrid,alpha+1))
      b = min(1d0,xg(igrid,beta)/xg(igrid,bound))
*
      fL = 0d0
      if(alpha.eq.beta) fL = 1d0
*
*     Initialize Integrals
*
      do inf=3,6
         do k=1,4
            do wipt=0,ipt
               SC2zm(igrid,inf,k,wipt,beta,alpha) = 0d0
               SCLzm(igrid,inf,k,wipt,beta,alpha) = 0d0
               SC3zm(igrid,inf,k,wipt,beta,alpha) = 0d0
            enddo
         enddo
*
*     Variables needed for wrapping the integrand functions
*
         wnf    = inf
         walpha = alpha
         wbeta  = beta
*
*     Precompute integrals
*
         if(ipt.ge.1)then
            wipt = 1
*
            sf = 1
            k = 1
            C2g1R = dgauss(integrandsSIAzm,a,b,eps(wipt))
            k = 3
            C2ns1RS = dgauss(integrandsSIAzm,a,b,eps(wipt))
            C2ns1L  = C2NS1TC(a)
*
            sf = 2
            k = 1
            CLg1R = dgauss(integrandsSIAzm,a,b,eps(wipt))
            k = 3
            CLns1RS = dgauss(integrandsSIAzm,a,b,eps(wipt))
*
            sf = 3
            k = 3
            C3ns1RS = dgauss(integrandsSIAzm,a,b,eps(wipt))
            C3ns1L  = C3NS1TC(a)
         endif
c$$$         if(ipt.ge.2)then
c$$$            wipt = 2
c$$$*
c$$$            sf = 1
c$$$            k = 1
c$$$            C2g2R    = dgauss(integrandsSIAzm,a,b,eps(wipt))
c$$$            C2g2L    = C2G2TC(a,1)
c$$$            k = 2
c$$$            C2ps2R   = dgauss(integrandsSIAzm,a,b,eps(wipt))
c$$$            k = 3
c$$$            C2nsp2RS = dgauss(integrandsSIAzm,a,b,eps(wipt))
c$$$            C2nsp2L  = C2NSP2TC(a,inf)
c$$$            k = 4
c$$$            C2nsm2RS = dgauss(integrandsSIAzm,a,b,eps(wipt))
c$$$            C2nsm2L  = C2NSM2TC(a,inf)
c$$$*
c$$$            sf = 2
c$$$            k = 1
c$$$            CLg2R    = dgauss(integrandsSIAzm,a,b,eps(wipt))
c$$$            k = 2
c$$$            CLps2R   = dgauss(integrandsSIAzm,a,b,eps(wipt))
c$$$            k = 3
c$$$            CLnsp2RS = dgauss(integrandsSIAzm,a,b,eps(wipt))
c$$$            CLnsp2L  = CLNSP2TC(a)
c$$$            k = 4
c$$$            CLnsm2RS = dgauss(integrandsSIAzm,a,b,eps(wipt))
c$$$            CLnsm2L  = CLNSM2TC(a)
c$$$*
c$$$            sf = 3
c$$$            k = 2
c$$$            C3ps2R   = dgauss(integrandsSIAzm,a,b,eps(wipt))
c$$$            k = 3
c$$$            C3nsp2RS = dgauss(integrandsSIAzm,a,b,eps(wipt))
c$$$            C3nsp2L  = C3NSP2TC(a,inf)
c$$$            k = 4
c$$$            C3nsm2RS = dgauss(integrandsSIAzm,a,b,eps(wipt))
c$$$            C3nsm2L  = C3NSM2TC(a,inf)
c$$$         endif
*
         do k=1,4
*
*     LO
*
*     C2
            C2L(k,0) = 0d0
            if(k.eq.3.or.k.eq.4) C2L(k,0) = 1d0
            integC2(0) = 0d0
*     CL
            CLL(k,0) = 0d0
            integCL(0) = 0d0
*     C3
            C3L(k,0) = 0d0
            if(k.eq.3.or.k.eq.4) C3L(k,0) = 1d0
            integC3(0) = 0d0
*
*     NLO
*
            if(ipt.ge.1)then
*     Gluon
               if(k.eq.1)then
*     C2
                  C2L(k,1)   = 0d0
                  integC2(1) = C2g1R
*     CL
                  CLL(k,1)   = 0d0
                  integCL(1) = CLg1R
*     C3
                  C3L(k,1)   = 0d0
                  integC3(1) = 0d0
*     Pure-Singlet
               elseif(k.eq.2)then
*     C2
                  C2L(k,1)   = 0d0
                  integC2(1) = 0d0
*     CL
                  CLL(k,1)   = 0d0
                  integCL(1) = 0d0
*     C3
                  C3L(k,1)   = 0d0
                  integC3(1) = 0d0
*     Non-singlet-plus/minus
               elseif(k.eq.3.or.k.eq.4)then
*     C2
                  C2L(k,1)   = C2ns1L
                  integC2(1) = C2ns1RS
*     CL
                  CLL(k,1)   = 0d0
                  integCL(1) = CLns1RS
*     C3
                  C3L(k,1)   = C3ns1L
                  integC3(1) = C3ns1RS
               endif
            endif
c$$$*
c$$$*     NNLO
c$$$*
c$$$            if(ipt.ge.2)then
c$$$*     Gluon
c$$$               if(k.eq.1)then
c$$$*     C2
c$$$                  C2L(k,2)   = C2g2L
c$$$                  integC2(2) = C2g2R
c$$$*     CL
c$$$                  CLL(k,2)   = 0d0
c$$$                  integCL(2) = CLg2R
c$$$*     C3
c$$$                  C3L(k,2)   = 0d0
c$$$                  integC3(2) = 0d0
c$$$*     Pure-Singlet
c$$$               elseif(k.eq.2)then
c$$$*     C2
c$$$                  C2L(k,2)   = 0d0
c$$$                  integC2(2) = C2ps2R
c$$$*     CL
c$$$                  CLL(k,2)   = 0d0
c$$$                  integCL(2) = CLps2R
c$$$*     C3
c$$$                  C3L(k,2)   = 0d0
c$$$                  integC3(2) = 0d0
c$$$*     Non-singlet-plus
c$$$               elseif(k.eq.3)then
c$$$*     C2
c$$$                  C2L(k,2)   = C2nsp2L
c$$$                  integC2(2) = C2nsp2RS
c$$$*     CL
c$$$                  CLL(k,2)   = CLnsp2L
c$$$                  integCL(2) = CLnsp2RS
c$$$*     C3
c$$$                  C3L(k,2)   = C3nsp2L
c$$$                  integC3(2) = C3nsp2RS
c$$$*     Non-singlet-minus
c$$$               elseif(k.eq.4)then
c$$$*     C2
c$$$                  C2L(k,2)   = C2nsm2L
c$$$                  integC2(2) = C2nsm2RS
c$$$*     CL
c$$$                  CLL(k,2)   = CLnsm2L
c$$$                  integCL(2) = CLnsm2RS
c$$$*     C3
c$$$                  C3L(k,2)   = C3nsm2L
c$$$                  integC3(2) = C3nsm2RS
c$$$               endif
c$$$            endif
*
*     Integrals
*
            do wipt=0,ipt
               SC2zm(igrid,inf,k,wipt,beta,alpha) = integC2(wipt) 
     1                                            + C2L(k,wipt) * fL
               SCLzm(igrid,inf,k,wipt,beta,alpha) = integCL(wipt) 
     1                                            + CLL(k,wipt) * fL
               SC3zm(igrid,inf,k,wipt,beta,alpha) = integC3(wipt) 
     1                                            + C3L(k,wipt) * fL
            enddo
         enddo
      enddo
*
      return
      end
