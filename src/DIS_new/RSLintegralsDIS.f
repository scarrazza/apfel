************************************************************************
*
*     RSLintegralsDIS.f:
*
*     This routine evaluates and dump on the common integrals.h once and 
*     for all for the pair (alpha,beta) the integral of the sum of Regular
*     and Singular part plus the Local part of all the coefficient functions
*     for all the needed orders in QCD.
*     The index sf runs like that:
*
*     sf  =   1   2   3
*             C2  CL  C3
*
*     while the index kk runs like that:
*
*     kk  =    1            2            3
*            gluon     pure-singlet  non-singlet
*
************************************************************************
      subroutine RSLintegralsDIS(beta,alpha)
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/wrapDIS.h"
      include "../commons/grid.h"
      include "../commons/integralsDIS.h"
**
*     Input Variables
*
      integer beta,alpha
**
*     Internal Variables
*
      integer bound,inf
      double precision C2L(3,0:2),CLL(3,0:2),C3L(3,0:2)
      double precision fL
      double precision dgauss,a,b,eps(2)
      double precision integrandsDIS
      double precision integC2(0:2),integCL(0:2),integC3(0:2)
      double precision C2ns1RS,C2ns1L,C2g1R,CLns1RS,CLg1R,C3ns1RS,C3ns1L
      double precision C2g2R,C2g2L,C2ps2R,C2ns2RS,C2ns2L,CLg2R,CLns2RS
      double precision CLps2R,CLns2L,C3ps2R,C3ns2RS,C3ns2L
      double precision C2NS1C,C3NS1C
      double precision C2G2C,C2NSP2C,CLNSP2C,C3NSP2C
      external integrandsDIS
c      data eps / 5d-8, 1d-3 /
      data eps / 5d-8, 1d-5 /
c      data eps / 1d-7, 1d-5 /
*
*     Initialize Integrals
*
      do inf=3,6
         do k=1,3
            do wipt=0,ipt
               SC2(igrid,inf,k,wipt,beta,alpha) = 0d0
               SCL(igrid,inf,k,wipt,beta,alpha) = 0d0
               SC3(igrid,inf,k,wipt,beta,alpha) = 0d0
            enddo
         enddo

*
*     Adjustment of the bounds of the integrals
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
            sf = 1
            k = 1
            C2g1R = dgauss(integrandsDIS,a,b,eps(wipt))
            k = 3
            C2ns1RS = dgauss(integrandsDIS,a,b,eps(wipt))
            C2ns1L  = C2NS1C(a)
            sf = 2
*
            k = 1
            CLg1R = dgauss(integrandsDIS,a,b,eps(wipt))
            k = 3
            CLns1RS = dgauss(integrandsDIS,a,b,eps(wipt))
            sf = 3
*
            k = 3
            C3ns1RS = dgauss(integrandsDIS,a,b,eps(wipt))
            C3ns1L  = C3NS1C(a)
         endif
         if(ipt.ge.2)then
            wipt = 2
            sf = 1
            k = 1
            C2g2R   = dgauss(integrandsDIS,a,b,eps(wipt))
            C2g2L   = C2G2C(a,1)
            k = 2
            C2ps2R  = dgauss(integrandsDIS,a,b,eps(wipt))
            k = 3
            C2ns2RS = dgauss(integrandsDIS,a,b,eps(wipt))
            C2ns2L  = C2NSP2C(a,inf)
*
            sf = 2
            k = 1
            CLg2R   = dgauss(integrandsDIS,a,b,eps(wipt))
            k = 2
            CLps2R  = dgauss(integrandsDIS,a,b,eps(wipt))
            k = 3
            CLns2RS = dgauss(integrandsDIS,a,b,eps(wipt))
            CLns2L  = CLNSP2C(a)
*
            sf = 3
            k = 2
            C3ps2R  = dgauss(integrandsDIS,a,b,eps(wipt))
            k = 3
            C3ns2RS = dgauss(integrandsDIS,a,b,eps(wipt))
            C3ns2L  = C3NSP2C(a,inf)
         endif
*
         do k=1,3
*
*     LO
*
*     C2
            C2L(k,0) = 0d0
            if(k.eq.3) C2L(k,0) = 1d0
            integC2(0) = 0d0
*     CL
            CLL(k,0) = 0d0
            integCL(0) = 0d0
*     C3
            C3L(k,0) = 0d0
            if(k.eq.3) C3L(k,0) = 1d0
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
*     Non-singlet
            elseif(k.eq.3)then
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
*
*     NNLO
*
         if(ipt.ge.2)then
*     Gluon
            if(k.eq.1)then
*     C2
               C2L(k,2)   = C2g2L
               integC2(2) = C2g2R
*     CL
               CLL(k,2)   = 0d0
               integCL(2) = CLg2R
*     C3
               C3L(k,2)   = 0d0
               integC3(2) = 0d0
*     Pure-Singlet
            elseif(k.eq.2)then
*     C2
               C2L(k,2)   = 0d0
               integC2(2) = C2ps2R
*     CL
               CLL(k,2)   = 0d0
               integCL(2) = CLps2R
*     C3
               C3L(k,2)   = 0d0
               integC3(2) = C2ps2R
*     Non-singlet
            elseif(k.eq.3)then
*     C2
               C2L(k,2)   = C2ns2L
               integC2(2) = C2ns2RS
*     CL
               CLL(k,2)   = CLns2L
               integCL(2) = CLns2RS
*     C3
               C3L(k,2)   = C3ns2L
               integC3(2) = C3ns2RS
            endif
         endif



c$$$*
c$$$*     NNLO
c$$$*
c$$$         if(ipt.ge.2)then
c$$$*     Plus
c$$$            if(k.eq.1)then
c$$$               PL(k,2)  = ns2Lp
c$$$               integ(2) = ns2RSp
c$$$*     Minus
c$$$            elseif(k.eq.2)then
c$$$               PL(k,2) = ns2Lm
c$$$               integ(2) = ns2RSm
c$$$*     Valence
c$$$            elseif(k.eq.3)then
c$$$               PL(k,2) = ns2Lm
c$$$               integ(2) = ns2RSv
c$$$*     Quark-Quark
c$$$            elseif(k.eq.4)then
c$$$               PL(k,2) = ns2Lp
c$$$               integ(2) = qq2RS
c$$$*     Quark-Gluon
c$$$            elseif(k.eq.5)then
c$$$               PL(k,2) = 0d0
c$$$               integ(2) = qg2RS
c$$$*     Gluon-Quark
c$$$            elseif(k.eq.6)then
c$$$               PL(k,2) = 0d0
c$$$               integ(2) = gq2RS
c$$$*     Gluon-Gluon
c$$$            elseif(k.eq.7)then
c$$$               PL(k,2) = gg2L
c$$$               integ(2) = gg2RS
c$$$            endif
c$$$         endif
*
*     Integrals
*
            do wipt=0,ipt
               SC2(igrid,inf,k,wipt,beta,alpha) = integC2(wipt) 
     1                                          + C2L(k,wipt) * fL
               SCL(igrid,inf,k,wipt,beta,alpha) = integCL(wipt) 
     1                                          + CLL(k,wipt) * fL
               SC3(igrid,inf,k,wipt,beta,alpha) = integC3(wipt) 
     1                                          + C3L(k,wipt) * fL
            enddo
c$$$*
c$$$*     In case of muR.ne.muF...
c$$$*     Eq. (2.8) of hep-ph/0408244
c$$$*
c$$$         if(kren.ne.1d0)then
c$$$            lnk = - dlog(kren)
c$$$            if(ipt.eq.1)then
c$$$               SP(igrid,nf,k,1,beta,alpha) = SP(igrid,nf,k,1,beta,alpha) 
c$$$     1         - beta0(nf) * lnk * SP(igrid,nf,k,0,beta,alpha) 
c$$$            elseif(ipt.eq.2)then
c$$$               SP(igrid,nf,k,2,beta,alpha) = SP(igrid,nf,k,2,beta,alpha) 
c$$$     1         - 2d0 * beta0(nf) * lnk * SP(igrid,nf,k,1,beta,alpha)
c$$$     2         - ( beta1(nf) - beta0(nf)**2d0 * lnk ) 
c$$$     3         * lnk * SP(igrid,nf,k,0,beta,alpha)
c$$$               SP(igrid,nf,k,1,beta,alpha) = SP(igrid,nf,k,1,beta,alpha) 
c$$$     1         - beta0(nf) * lnk * SP(igrid,nf,k,0,beta,alpha) 
c$$$            endif
c$$$         endif
         enddo
      enddo
*
      return
      end
