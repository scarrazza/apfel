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
      include "../commons/coeffhqmellin.h"
      include "../commons/integralsDIS.h"
      include "../commons/MassScheme.h"
      include "../commons/Nf_FF.h"
**
*     Input Variables
*
      integer beta,alpha
**
*     Internal Variables
*
      integer bound,inf,ixi
      double precision C2L(3,0:2),CLL(3,0:2),C3L(3,0:2)
      double precision fL
      double precision dgauss,a,b,eps(2)
      double precision integrandsDISzm,integrandsDISm,integrandsDISm0
      double precision cm22q_adler
      double precision integC2(0:2),integCL(0:2),integC3(0:2)
      double precision C2ns1RS,C2ns1L,C2g1R,CLns1RS,CLg1R,C3ns1RS,C3ns1L
      double precision C2g2R,C2g2L,C2ps2R,C2ns2RS,C2ns2L,CLg2R,CLns2RS
      double precision CLps2R,CLns2L,C3ps2R,C3ns2RS,C3ns2L
      double precision C2g1R0,C2g1RQ,CLg1R0
      double precision C2g2R0,C2g2RQ,C2g2RQ2,C2g2RF,C2g2RQF
      double precision C2ps2R0,C2ps2RQ,C2ps2RQ2,C2ps2RF,C2ps2RQF
      double precision C2ns2RS0,C2ns2RSQ,C2ns2RSQ2
      double precision C2ns2L0,C2ns2LQ,C2ns2LQ2
      double precision CLg2R0,CLg2RQ,CLg2RF
      double precision CLps2R0,CLps2RQ,CLps2RF
      double precision CLns2R0,CLns2RQ
      double precision C2NS1C,C3NS1C
      double precision C2G2C,C2NSP2C,CLNSP2C,C3NSP2C
      double precision C2NS2CM0_A0,C2NS2CM0_AQ,C2NS2CM0_AQ2
      double precision kQF2,lnQ,lnF,lnQ2,lnQF
      external integrandsDISzm,integrandsDISm,integrandsDISm0
      external cm22q_adler
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
*
*     ZM-VFNS
*
      if(MassScheme.eq."ZM-VFNS".or.MassScheme(1:5).eq."FONLL".or.
     1   MassScheme(1:4).eq."FFNS".or.MassScheme(1:4).eq."FFN0")then
         do inf=3,6
            do k=1,3
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
               C2g1R = dgauss(integrandsDISzm,a,b,eps(wipt))
               k = 3
               C2ns1RS = dgauss(integrandsDISzm,a,b,eps(wipt))
               C2ns1L  = C2NS1C(a)
*
               sf = 2
               k = 1
               CLg1R = dgauss(integrandsDISzm,a,b,eps(wipt))
               k = 3
               CLns1RS = dgauss(integrandsDISzm,a,b,eps(wipt))
*
               sf = 3
               k = 3
               C3ns1RS = dgauss(integrandsDISzm,a,b,eps(wipt))
               C3ns1L  = C3NS1C(a)
            endif
            if(ipt.ge.2)then
               wipt = 2
*
               sf = 1
               k = 1
               C2g2R   = dgauss(integrandsDISzm,a,b,eps(wipt))
               C2g2L   = C2G2C(a,1)
               k = 2
               C2ps2R  = dgauss(integrandsDISzm,a,b,eps(wipt))
               k = 3
               C2ns2RS = dgauss(integrandsDISzm,a,b,eps(wipt))
               C2ns2L  = C2NSP2C(a,inf)
*
               sf = 2
               k = 1
               CLg2R   = dgauss(integrandsDISzm,a,b,eps(wipt))
               k = 2
               CLps2R  = dgauss(integrandsDISzm,a,b,eps(wipt))
               k = 3
               CLns2RS = dgauss(integrandsDISzm,a,b,eps(wipt))
               CLns2L  = CLNSP2C(a)
*
               sf = 3
               k = 2
               C3ps2R  = dgauss(integrandsDISzm,a,b,eps(wipt))
               k = 3
               C3ns2RS = dgauss(integrandsDISzm,a,b,eps(wipt))
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
                     integC3(2) = 0d0
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
*
*     Integrals
*
               do wipt=0,ipt
                  SC2zm(igrid,inf,k,wipt,beta,alpha) = integC2(wipt) 
     1                                               + C2L(k,wipt) * fL
                  SCLzm(igrid,inf,k,wipt,beta,alpha) = integCL(wipt) 
     1                                               + CLL(k,wipt) * fL
                  SC3zm(igrid,inf,k,wipt,beta,alpha) = integC3(wipt) 
     1                                               + C3L(k,wipt) * fL
               enddo
            enddo
         enddo
      endif
*
*     FFNS
*
      if(MassScheme(1:4).eq."FFNS".or.MassScheme(1:5).eq."FONLL")then
         do ixi=1,nxi
            do k=1,3
               do wipt=0,ipt
                  SC2m(igrid,ixi,k,wipt,beta,alpha) = 0d0
                  SCLm(igrid,ixi,k,wipt,beta,alpha) = 0d0
                  SC3m(igrid,ixi,k,wipt,beta,alpha) = 0d0
               enddo
            enddo
*
*     Variables needed for wrapping the integrand functions
*
            wnf    = Nf_FF
            wixi   = ixi
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
               C2g1R = dgauss(integrandsDISm,a,b,eps(wipt))
*
               sf = 2
               k = 1
               CLg1R = dgauss(integrandsDISm,a,b,eps(wipt))
*
               sf = 3
               k = 3
               C3ns1RS = dgauss(integrandsDISzm,a,b,eps(wipt))
               C3ns1L  = C3NS1C(a)
            endif
            if(ipt.ge.2)then
               wipt = 2
*
               sf = 1
               k = 1
               C2g2R   = dgauss(integrandsDISm,a,b,eps(wipt))
               k = 2
               C2ps2R  = dgauss(integrandsDISm,a,b,eps(wipt))
               k = 3
               C2ns2RS = dgauss(integrandsDISm,a,b,eps(wipt))
               C2ns2L  = - dgauss(cm22q_adler,0d0,1d0,eps(wipt))
*
               sf = 2
               k = 1
               CLg2R   = dgauss(integrandsDISm,a,b,eps(wipt))
               k = 2
               CLps2R  = dgauss(integrandsDISm,a,b,eps(wipt))
               k = 3
               CLns2RS = dgauss(integrandsDISm,a,b,eps(wipt))
*
               sf = 3
               k = 2
               C3ps2R  = dgauss(integrandsDISm,a,b,eps(wipt))
               k = 3
               C3ns2RS = dgauss(integrandsDISm,a,b,eps(wipt))
               C3ns2L  = C3NSP2C(a,inf)
            endif
*
            do k=1,3
*
*     LO
*
*     C2
               C2L(k,0) = 0d0
               integC2(0) = 0d0
*     CL
               integCL(0) = 0d0
*     C3
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
                     C2L(k,1) = 0d0
                     integC2(1) = C2g1R
*     CL
                     integCL(1) = CLg1R
*     C3
                     C3L(k,1)   = 0d0
                     integC3(1) = 0d0
*     Pure-Singlet
                  elseif(k.eq.2)then
*     C2
                     C2L(k,1) = 0d0
                     integC2(1) = 0d0
*     CL
                     integCL(1) = 0d0
*     C3
                     C3L(k,1)   = 0d0
                     integC3(1) = 0d0
*     Non-singlet
                  elseif(k.eq.3)then
*     C2
                     C2L(k,1) = 0d0
                     integC2(1) = 0d0
*     CL
                     integCL(1) = 0d0
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
                     C2L(k,2) = 0d0
                     integC2(2) = C2g2R
*     CL
                     integCL(2) = CLg2R
*     C3
                     C3L(k,2)   = 0d0
                     integC3(2) = 0d0
*     Pure-Singlet
                  elseif(k.eq.2)then
*     C2
                     C2L(k,2) = 0d0
                     integC2(2) = C2ps2R
*     CL
                     integCL(2) = CLps2R
*     C3
                     C3L(k,2)   = 0d0
                     integC3(2) = 0d0
*     Non-singlet
                  elseif(k.eq.3)then
*     C2
                     C2L(k,2) = C2ns2L
                     integC2(2) = C2ns2RS
*     CL
                     integCL(2) = CLns2RS
*     C3
                     C3L(k,2)   = C3ns2L
                     integC3(2) = C3ns2RS
                  endif
               endif
*
*     Integrals
*

               do wipt=0,ipt
                  SC2m(igrid,ixi,k,wipt,beta,alpha) = integC2(wipt) 
     1                                              + C2L(k,wipt) * fL
                  SCLm(igrid,ixi,k,wipt,beta,alpha) = integCL(wipt) 
                  SC3m(igrid,ixi,k,wipt,beta,alpha) = integC3(wipt)
     1                                              + C3L(k,wipt) * fL
               enddo
            enddo
         enddo
      endif
*
*     Massive zero FFNS needed for the FONLL scheme
*
      if(MassScheme(1:4).eq."FFN0".or.MassScheme(1:5).eq."FONLL")then
*
*     Variables needed for wrapping the integrand functions
*
         wnf    = Nf_FF
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
            wl = 1
            C2g1R0 = dgauss(integrandsDISm0,a,b,eps(wipt))
            wl = 2
            C2g1RQ = dgauss(integrandsDISm0,a,b,eps(wipt))
*
            sf = 2
            k = 1
            wl = 1
            CLg1R0 = dgauss(integrandsDISm0,a,b,eps(wipt))
*
            sf = 3
            k = 3
            C3ns1RS = dgauss(integrandsDISzm,a,b,eps(wipt))
            C3ns1L  = C3NS1C(a)
         endif
         if(ipt.ge.2)then
            wipt = 2
*
            sf = 1
            k = 1
            wl = 1
            C2g2R0  = dgauss(integrandsDISm0,a,b,eps(wipt))
            wl = 2
            C2g2RQ  = dgauss(integrandsDISm0,a,b,eps(wipt))
            wl = 3
            C2g2RQ2 = dgauss(integrandsDISm0,a,b,eps(wipt))
            wl = 4
            C2g2RF  = dgauss(integrandsDISm0,a,b,eps(wipt))
            wl = 5
            C2g2RQF = dgauss(integrandsDISm0,a,b,eps(wipt))
*
            k = 2
            wl = 1
            C2ps2R0  = dgauss(integrandsDISm0,a,b,eps(wipt))
            wl = 2
            C2ps2RQ  = dgauss(integrandsDISm0,a,b,eps(wipt))
            wl = 3
            C2ps2RQ2 = dgauss(integrandsDISm0,a,b,eps(wipt))
            wl = 4
            C2ps2RF  = dgauss(integrandsDISm0,a,b,eps(wipt))
            wl = 5
            C2ps2RQF = dgauss(integrandsDISm0,a,b,eps(wipt))
*
            k = 3
            wl = 1
            C2ns2RS0  = dgauss(integrandsDISm0,a,b,eps(wipt))
            C2ns2L0   = C2NS2CM0_A0(a)
            wl = 2
            C2ns2RSQ  = dgauss(integrandsDISm0,a,b,eps(wipt))
            C2ns2LQ   = C2NS2CM0_AQ(a)
            wl = 3
            C2ns2RSQ2 = dgauss(integrandsDISm0,a,b,eps(wipt))
            C2ns2LQ2  = C2NS2CM0_AQ2(a)
*
            sf = 2
            k = 1
            wl = 1
            CLg2R0 = dgauss(integrandsDISm0,a,b,eps(wipt))
            wl = 2
            CLg2RQ = dgauss(integrandsDISm0,a,b,eps(wipt))
            wl = 4
            CLg2RF = dgauss(integrandsDISm0,a,b,eps(wipt))
*
            k = 2
            wl = 1
            CLps2R0 = dgauss(integrandsDISm0,a,b,eps(wipt))
            wl = 2
            CLps2RQ = dgauss(integrandsDISm0,a,b,eps(wipt))
            wl = 4
            CLps2RF = dgauss(integrandsDISm0,a,b,eps(wipt))
*
            k = 3
            wl = 1
            CLns2R0 = dgauss(integrandsDISm0,a,b,eps(wipt))
            wl = 2
            CLns2RQ = dgauss(integrandsDISm0,a,b,eps(wipt))
*
            sf = 3
            k = 2
            C3ps2R  = dgauss(integrandsDISm0,a,b,eps(wipt))
            k = 3
            C3ns2RS = dgauss(integrandsDISm0,a,b,eps(wipt))
            C3ns2L  = C3NSP2C(a,inf)
         endif
*
         do ixi=1,nxi
            kQF2 = 1d0 ! Q2 / muF2
            lnQ  = dlog(xigrid(ixi))
            lnF  = dlog(kQF2) - lnQ
            lnQ2 = lnQ * lnQ
            lnQF = lnQ * lnF
            do k=1,3
               do wipt=0,ipt
                  SC2m0(igrid,ixi,k,wipt,beta,alpha) = 0d0
                  SCLm0(igrid,ixi,k,wipt,beta,alpha) = 0d0
                  SC3m0(igrid,ixi,k,wipt,beta,alpha) = 0d0
               enddo
            enddo
*
            do k=1,3
*
*     LO
*
*     C2
               integC2(0) = 0d0
*     CL
               integCL(0) = 0d0
*     C3
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
                     C2L(k,1) = 0d0
                     integC2(1) = C2g1R0 + C2g1RQ * lnQ
*     CL
                     integCL(1) = CLg1R0
*     C3
                     C3L(k,1)   = 0d0
                     integC3(1) = 0d0
*     Pure-Singlet
                  elseif(k.eq.2)then
*     C2
                     C2L(k,1) = 0d0
                     integC2(1) = 0d0
*     CL
                     integCL(1) = 0d0
*     C3
                     C3L(k,1)   = 0d0
                     integC3(1) = 0d0
*     Non-singlet
                  elseif(k.eq.3)then
*     C2
                     C2L(k,1) = 0d0
                     integC2(1) = 0d0
*     CL
                     integCL(1) = 0d0
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
                     C2L(k,2) = 0d0
                     integC2(2) = C2g2R0 + C2g2RQ * lnQ + C2g2RQ2 * lnQ2 
     1                          + C2g2RF * lnF + C2g2RQF * lnQF
*     CL
                     integCL(2) = CLg2R0 + CLg2RQ * lnQ + CLg2RF * lnF
*     C3
                     C3L(k,2)   = 0d0
                     integC3(2) = 0d0
*     Pure-Singlet
                  elseif(k.eq.2)then
*     C2
                     C2L(k,2) = 0d0
                     integC2(2) = C2ps2R0 + C2ps2RQ * lnQ 
     1                          + C2ps2RQ2 * lnQ2 + C2ps2RF * lnF
     2                          + C2ps2RQF * lnQF
*     CL
                     integCL(2) = CLps2R0 + CLps2RQ * lnQ 
     1                          + CLps2RF * lnF
*     C3
                     C3L(k,2)   = 0d0
                     integC3(2) = 0d0
*     Non-singlet
                  elseif(k.eq.3)then
*     C2
                     C2L(k,2) = C2ns2L0 + C2ns2LQ * lnQ 
     1                        + C2ns2LQ2 * lnQ2
                     integC2(2) = C2ns2RS0 + C2ns2RSQ * lnQ 
     1                          + C2ns2RSQ2 * lnQ2
*     CL
                     integCL(2) = CLns2R0 + CLns2RQ * lnQ
*     C3
                     C3L(k,2)   = C3ns2L
                     integC3(2) = C3ns2RS
                  endif
               endif
*
*     Integrals
*
               do wipt=0,ipt
                  SC2m0(igrid,ixi,k,wipt,beta,alpha) = integC2(wipt) 
     1                                               + C2L(k,wipt) * fL
                  SCLm0(igrid,ixi,k,wipt,beta,alpha) = integCL(wipt) 
                  SC3m0(igrid,ixi,k,wipt,beta,alpha) = integC3(wipt)
     1                                               + C3L(k,wipt) * fL
               enddo
            enddo
         enddo
      endif
*
      return
      end
