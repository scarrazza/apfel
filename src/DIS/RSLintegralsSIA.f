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
c      include "../commons/ipt.h"
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
      integer bound,inf,ipt_ZM
      double precision C2L(0:2),CLL(0:2),C3L(0:2)
      double precision fL
      double precision dgauss,a,b,eps(2)
      double precision integrandsSIAzm
      double precision integC2(0:2),integCL(0:2),integC3(0:2)
      double precision C2ns1RS,C2ns1L,C2g1R,CLns1RS,CLg1R,C3ns1RS,C3ns1L
      double precision C2g2R,CLg2R,C2ps2R,CLps2R,C3ps2R
      double precision C2nsp2RS,CLnsp2RS,C3nsp2RS
      double precision C2nsp2L,C3nsp2L
      double precision C2nsm2RS,CLnsm2RS,C3nsm2RS
      double precision C2nsm2L,C3nsm2L
      double precision C2NS1TC,C3NS1TC
      double precision C2NSP2TC,C3NSP2TC
      external integrandsSIAzm
c      data eps / 5d-8, 1d-3 /
c      data eps / 5d-8, 1d-5 /
      data eps / 1d-5, 1d-3 /
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
*     The ZM scheme is cheap enough to compute all orders anyway. This
*     is useful because in some cases one may need different absolute
*     orderes at the same time, like with F2 and FL.
*
      ipt_ZM = 2
*
*     Initialize Integrals
*
      do inf=3,6
         do k=1,4
            do wipt=0,ipt_ZM
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
         if(ipt_ZM.ge.1)then
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
         if(ipt_ZM.ge.2)then
            wipt = 2
*
            sf = 1
            k = 1
            C2g2R    = dgauss(integrandsSIAzm,a,b,eps(wipt))
            k = 2
            C2ps2R   = dgauss(integrandsSIAzm,a,b,eps(wipt))
            k = 3
            C2nsp2RS = dgauss(integrandsSIAzm,a,b,eps(wipt))
            C2nsp2L  = C2NSP2TC(a,inf)
            k = 4
            C2nsm2RS = C2nsp2RS
            C2nsm2L  = C2nsp2L
*
            sf = 2
            k = 1
            CLg2R    = dgauss(integrandsSIAzm,a,b,eps(wipt))
            k = 2
            CLps2R   = dgauss(integrandsSIAzm,a,b,eps(wipt))
            k = 3
            CLnsp2RS = dgauss(integrandsSIAzm,a,b,eps(wipt))
            k = 4
            CLnsm2RS = CLnsp2RS
*
            sf = 3
            k = 2
            C3ps2R   = dgauss(integrandsSIAzm,a,b,eps(wipt))
            k = 3
            C3nsp2RS = dgauss(integrandsSIAzm,a,b,eps(wipt))
            C3nsp2L  = C3NSP2TC(a,inf)
            k = 4
            C3nsm2RS = C3nsp2RS
            C3nsm2L  = C3nsp2L
         endif
*
         do k=1,4
*
*     LO
*
*     C2
            C2L(0) = 0d0
            if(k.eq.3.or.k.eq.4) C2L(0) = 1d0
            integC2(0) = 0d0
*     CL
            CLL(0) = 0d0
            integCL(0) = 0d0
*     C3
            C3L(0) = 0d0
            if(k.eq.3.or.k.eq.4) C3L(0) = 1d0
            integC3(0) = 0d0
*
*     NLO
*
            if(ipt_ZM.ge.1)then
*     Gluon
               if(k.eq.1)then
*     C2
                  C2L(1)   = 0d0
                  integC2(1) = C2g1R
*     CL
                  CLL(1)   = 0d0
                  integCL(1) = CLg1R
*     C3
                  C3L(1)   = 0d0
                  integC3(1) = 0d0
*     Pure-Singlet
               elseif(k.eq.2)then
*     C2
                  C2L(1)   = 0d0
                  integC2(1) = 0d0
*     CL
                  CLL(1)   = 0d0
                  integCL(1) = 0d0
*     C3
                  C3L(1)   = 0d0
                  integC3(1) = 0d0
*     Non-singlet-plus/minus
               elseif(k.eq.3.or.k.eq.4)then
*     C2
                  C2L(1)   = C2ns1L
                  integC2(1) = C2ns1RS
*     CL
                  CLL(1)   = 0d0
                  integCL(1) = CLns1RS
*     C3
                  C3L(1)   = C3ns1L
                  integC3(1) = C3ns1RS
               endif
            endif
*
*     NNLO
*
            if(ipt_ZM.ge.2)then
*     Gluon
               if(k.eq.1)then
*     C2
                  C2L(2)   = 0d0
                  integC2(2) = C2g2R
*     CL
                  CLL(2)   = 0d0
                  integCL(2) = CLg2R
*     C3
                  C3L(2)   = 0d0
                  integC3(2) = 0d0
*     Pure-Singlet
               elseif(k.eq.2)then
*     C2
                  C2L(2)   = 0d0
                  integC2(2) = C2ps2R
*     CL
                  CLL(2)   = 0d0
                  integCL(2) = CLps2R
*     C3
                  C3L(2)   = 0d0
                  integC3(2) = 0d0
*     Non-singlet-plus
               elseif(k.eq.3)then
*     C2
                  C2L(2)   = C2nsp2L
                  integC2(2) = C2nsp2RS
*     CL
                  CLL(2)   = 0d0
                  integCL(2) = CLnsp2RS
*     C3
                  C3L(2)   = C3nsp2L
                  integC3(2) = C3nsp2RS
*     Non-singlet-minus
               elseif(k.eq.4)then
*     C2
                  C2L(2)   = C2nsm2L
                  integC2(2) = C2nsm2RS
*     CL
                  CLL(2)   = 0d0
                  integCL(2) = CLnsm2RS
*     C3
                  C3L(2)   = C3nsm2L
                  integC3(2) = C3nsm2RS
               endif
            endif
*
*     Integrals
*
            do wipt=0,ipt_ZM
               SC2zm(igrid,inf,k,wipt,beta,alpha) = integC2(wipt) 
     1                                            + C2L(wipt) * fL
               SCLzm(igrid,inf,k,wipt,beta,alpha) = integCL(wipt) 
     1                                            + CLL(wipt) * fL
               SC3zm(igrid,inf,k,wipt,beta,alpha) = integC3(wipt) 
     1                                            + C3L(wipt) * fL
            enddo
         enddo
      enddo
*
      return
      end
