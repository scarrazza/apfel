************************************************************************
*
*     RSLintegralsDIS.f:
*
*     This routine evaluates and dump on the common integralsDIS.h once and 
*     for all for the pair (alpha,beta) the integral of the sum of Regular
*     and Singular part plus the Local part of all the coefficient functions
*     for all the needed orders in DIS.
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
      include "../commons/mass_scheme.h"
      include "../commons/MassRunning.h"
      include "../commons/ColorFactors.h"
      include "../commons/kfacQ.h"
      include "../commons/krenQ.h"
**
*     Input Variables
*
      integer beta,alpha
**
*     Internal Variables
*
      integer bound,inf,ixi,iptmx,ipt_FF
      double precision C2L(0:2),CLL(0:2),C3L(0:2)
      double precision fL,fL_CCm,dfL_CCm,w_int,dw_int
      double precision dgauss,a,b,eps(2),b_CCm
      double precision integrandsDISzm
      double precision integrandsDISNCm,integrandsDISNCm0
      double precision integrandsDISCCm,integrandsDISCCm0
      double precision cm22q_adler
      double precision CGCOLM0,lnC2,lnC3,CGM0
      double precision integC2(0:2),integCL(0:2),integC3(0:2)
      double precision C2ns1RS,C2ns1L,C2g1R,CLns1RS,CLg1R,C3ns1RS,C3ns1L
      double precision CLns1L,C3g1R
      double precision C2g2R,C2g2L,CLg2R,C2ps2R,CLps2R,C3ps2R
      double precision C2ns2RS,CLns2RS,C3ns2RS
      double precision C2ns2L,C3ns2L
      double precision C2nsp2RS,CLnsp2RS,C3nsp2RS
      double precision C2nsp2L,CLnsp2L,C3nsp2L
      double precision C2nsm2RS,CLnsm2RS,C3nsm2RS
      double precision C2nsm2L,CLnsm2L,C3nsm2L
      double precision C2g1R0,C2g1RQ,CLg1R0
      double precision C2g2R0,C2g2RQ,C2g2RQ2,C2g2RF,C2g2RQF
      double precision C2ps2R0,C2ps2RQ,C2ps2RQ2,C2ps2RF,C2ps2RQF
      double precision C2ns2RS0,C2ns2RSQ,C2ns2RSQ2
      double precision C2ns2L0,C2ns2LQ,C2ns2LQ2
      double precision CLg2R0,CLg2RQ,CLg2RF
      double precision CLps2R0,CLps2RQ,CLps2RF
      double precision CLns2R0,CLns2RQ
      double precision C2NS1C,C3NS1C
      double precision C2G2C
      double precision C2NSP2C,CLNSP2C,C3NSP2C,C2NSM2C,CLNSM2C,C3NSM2C
      double precision C2NS2CM0_A0,C2NS2CM0_AQ,C2NS2CM0_AQ2
      double precision c2ns1ccc,clns1ccc,c3ns1ccc
      double precision kQF2,lnkQF2,lnQ,lnF,lnQ2,lnQF
      double precision lambda,Rf,Rfun
      double precision h1
      external integrandsDISzm
      external integrandsDISNCm,integrandsDISNCm0
      external integrandsDISCCm,integrandsDISCCm0
      external cm22q_adler
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
*     Ration between Scale and factorization scale squared (to be put in a common)
*
      kQF2 = 1d0 / kfacQ               ! Q2 / muF2
      lnkQF2 = dlog(kQF2)
*
*     Initialize Integrals
*
*     ZM-VFNS
*
      if(MassScheme.eq."ZM-VFNS".or.MassScheme(1:5).eq."FONLL".or.
     1   MassScheme(1:4).eq."FFNS".or.MassScheme(1:4).eq."FFN0")then
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
               C2g2R    = dgauss(integrandsDISzm,a,b,eps(wipt))
               C2g2L    = C2G2C(a,1)
               k = 2
               C2ps2R   = dgauss(integrandsDISzm,a,b,eps(wipt))
               k = 3
               C2nsp2RS = dgauss(integrandsDISzm,a,b,eps(wipt))
               C2nsp2L  = C2NSP2C(a,inf)
               k = 4
               C2nsm2RS = dgauss(integrandsDISzm,a,b,eps(wipt))
               C2nsm2L  = C2NSM2C(a,inf)
*
               sf = 2
               k = 1
               CLg2R    = dgauss(integrandsDISzm,a,b,eps(wipt))
               k = 2
               CLps2R   = dgauss(integrandsDISzm,a,b,eps(wipt))
               k = 3
               CLnsp2RS = dgauss(integrandsDISzm,a,b,eps(wipt))
               CLnsp2L  = CLNSP2C(a)
               k = 4
               CLnsm2RS = dgauss(integrandsDISzm,a,b,eps(wipt))
               CLnsm2L  = CLNSM2C(a)
*
               sf = 3
               k = 2
               C3ps2R   = dgauss(integrandsDISzm,a,b,eps(wipt))
               k = 3
               C3nsp2RS = dgauss(integrandsDISzm,a,b,eps(wipt))
               C3nsp2L  = C3NSP2C(a,inf)
               k = 4
               C3nsm2RS = dgauss(integrandsDISzm,a,b,eps(wipt))
               C3nsm2L  = C3NSM2C(a,inf)
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
               if(ipt.ge.1)then
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
               if(ipt.ge.2)then
*     Gluon
                  if(k.eq.1)then
*     C2
                     C2L(2)   = C2g2L
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
                     CLL(2)   = CLnsp2L
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
                     CLL(2)   = CLnsm2L
                     integCL(2) = CLnsm2RS
*     C3
                     C3L(2)   = C3nsm2L
                     integC3(2) = C3nsm2RS
                  endif
               endif
*
*     Integrals
*
               do wipt=0,ipt
                  SC2zm(igrid,inf,k,wipt,beta,alpha) = integC2(wipt) 
     1                                               + C2L(wipt) * fL
                  SCLzm(igrid,inf,k,wipt,beta,alpha) = integCL(wipt) 
     1                                               + CLL(wipt) * fL
                  SC3zm(igrid,inf,k,wipt,beta,alpha) = integC3(wipt) 
     1                                               + C3L(wipt) * fL
               enddo
            enddo
         enddo
      endif
*
*     FFNS
*
      if(MassScheme(1:4).eq."FFNS".or.MassScheme(1:5).eq."FONLL")then
*
*     Neutral Current
*
*     In case of FONLL-B
*
         ipt_FF = ipt
         if(MassScheme.eq."FONLL-B".and.ipt.ge.1) ipt_FF = 2
*
         do ixi=1,nxir
            do k=1,3
               do wipt=1,ipt_FF
                  SC2mNC(igrid,ixi,k,wipt,beta,alpha) = 0d0
                  SCLmNC(igrid,ixi,k,wipt,beta,alpha) = 0d0
               enddo
            enddo
*
*     Variables needed for wrapping the integrand functions
*
            wnf    = Nf_FF
            wixi   = ixi * xistep
            walpha = alpha
            wbeta  = beta
*
*     Precompute integrals
*
            if(ipt_FF.ge.1)then
               wipt = 1
*
               sf = 1
               k = 1
               C2g1R = dgauss(integrandsDISNCm,a,b,eps(wipt))
*
               sf = 2
               k = 1
               CLg1R = dgauss(integrandsDISNCm,a,b,eps(wipt))
            endif
            if(ipt_FF.ge.2)then
*
               wipt = 2
*
               sf = 1
               k = 1
               C2g2R   = dgauss(integrandsDISNCm,a,b,eps(wipt))
               k = 2
               C2ps2R  = dgauss(integrandsDISNCm,a,b,eps(wipt))
               k = 3
               C2ns2RS = dgauss(integrandsDISNCm,a,b,eps(wipt))
*
*     Coefficient to be added to the NNLO non-singlet coefficient function
*     of F2 to fullfil the Adler sum rule (can be done only the first time).
*
               if(alpha.eq.0.and.beta.eq.0.and.igrid.eq.1)
     1              adler_coef(ixi) = 
     2              - dgauss(cm22q_adler,0d0,1d0,eps(wipt))
               C2ns2L  = adler_coef(ixi)
*
               sf = 2
               k = 1
               CLg2R   = dgauss(integrandsDISNCm,a,b,eps(wipt))
               k = 2
               CLps2R  = dgauss(integrandsDISNCm,a,b,eps(wipt))
               k = 3
               CLns2RS = dgauss(integrandsDISNCm,a,b,eps(wipt))
            endif
*
            do k=1,3
*
*     NLO
*
               if(ipt_FF.ge.1)then
*     Gluon
                  if(k.eq.1)then
*     C2
                     C2L(1) = 0d0
                     integC2(1) = C2g1R
*     CL
                     integCL(1) = CLg1R
*     Pure-Singlet
                  elseif(k.eq.2)then
*     C2
                     C2L(1) = 0d0
                     integC2(1) = 0d0
*     CL
                     integCL(1) = 0d0
*     Non-singlet
                  elseif(k.eq.3)then
*     C2
                     C2L(1) = 0d0
                     integC2(1) = 0d0
*     CL
                     integCL(1) = 0d0
                  endif
               endif
*
*     NNLO
*
               if(ipt_FF.ge.2)then
*     Gluon
                  if(k.eq.1)then
*     C2
                     C2L(2) = 0d0
                     integC2(2) = C2g2R
*     CL
                     integCL(2) = CLg2R
*     Pure-Singlet
                  elseif(k.eq.2)then
*     C2
                     C2L(2) = 0d0
                     integC2(2) = C2ps2R
*     CL
                     integCL(2) = CLps2R
*     Non-singlet
                  elseif(k.eq.3)then
*     C2
                     C2L(2) = C2ns2L
                     integC2(2) = C2ns2RS
*     CL
                     integCL(2) = CLns2RS
                  endif
               endif
*
*     Integrals
*
               do wipt=1,ipt_FF
                  SC2mNC(igrid,ixi,k,wipt,beta,alpha) = integC2(wipt) 
     1                                                + C2L(wipt) * fL
                  SCLmNC(igrid,ixi,k,wipt,beta,alpha) = integCL(wipt) 
               enddo
            enddo
         enddo
*
*     Charged Current
*
         iptmx = ipt
         if(ipt.gt.1) iptmx = 1
*
         do ixi=1,nxir
            do k=1,3
               do wipt=0,iptmx
                  SC2mCC(igrid,ixi,k,wipt,beta,alpha) = 0d0
                  SCLmCC(igrid,ixi,k,wipt,beta,alpha) = 0d0
                  SC3mCC(igrid,ixi,k,wipt,beta,alpha) = 0d0
               enddo
            enddo
*
            lambda = xigrid(ixi*xistep) / ( 1d0 + xigrid(ixi*xistep) )
*
            if(a.ge.lambda) cycle
*
            Rf     = RFun(xigrid(ixi*xistep),a/lambda)
            fL_CCm = w_int(inter_degree(igrid),alpha,
     1                     xg(igrid,beta)/lambda)
            if(mass_scheme.eq."MSbar")
     1           dfL_CCm = xg(igrid,beta) * dw_int(1,
     2           inter_degree(igrid),alpha,xg(igrid,beta)/lambda)
     3           / lambda
            b_CCm = min(lambda,xg(igrid,beta)/xg(igrid,bound))
*
*     Variables needed for wrapping the integrand functions
*
            wnf    = Nf_FF
            wixi   = ixi * xistep
            walpha = alpha
            wbeta  = beta
*
*     Precompute integrals
*
            if(iptmx.ge.1)then
               wipt = 1
*
               sf = 1
               k = 1
               C2g1R   = dgauss(integrandsDISCCm,a,b_CCm,eps(wipt))
               k = 3
               C2ns1RS = dgauss(integrandsDISCCm,a,b_CCm,eps(wipt))
               C2ns1L  = c2ns1ccc(Rf,xigrid(ixi*xistep),a/lambda)
*
               sf = 2
               k = 1
               CLg1R   = dgauss(integrandsDISCCm,a,b_CCm,eps(wipt))
               k = 3
               CLns1RS = dgauss(integrandsDISCCm,a,b_CCm,eps(wipt))
               CLns1L  = clns1ccc(Rf,xigrid(ixi*xistep),a/lambda)
*
               sf = 3
               k = 1
               C3g1R   = dgauss(integrandsDISCCm,a,b_CCm,eps(wipt))
               k = 3
               C3ns1RS = dgauss(integrandsDISCCm,a,b_CCm,eps(wipt))
               C3ns1L  = c3ns1ccc(Rf,xigrid(ixi*xistep),a/lambda)
            endif
*
            do k=1,3
*
*     LO
*
*     C2
               C2L(0)   = 0d0
               if(k.eq.3) C2L(0) = 1d0
               integC2(0) = 0d0
*     CL
               CLL(0)   = 0d0
               if(k.eq.3) CLL(0) = 1d0 - lambda
               integCL(0) = 0d0
*     C3
               C3L(0) = 0d0
               if(k.eq.3) C3L(0) = 1d0
               integC3(0) = 0d0
*
*     NLO
*
               if(iptmx.ge.1)then
*     Gluon
                  if(k.eq.1)then
*     C2
                     C2L(1) = 0d0
                     integC2(1) = C2g1R
*     CL
                     CLL(1) = 0d0
                     integCL(1) = CLg1R
*     C3
                     C3L(1) = 0d0
                     integC3(1) = C3g1R
*     Pure-Singlet
                  elseif(k.eq.2)then
*     C2
                     C2L(1) = 0d0
                     integC2(1) = 0d0
*     CL
                     CLL(1) = 0d0
                     integCL(1) = 0d0
*     C3
                     C3L(1) = 0d0
                     integC3(1) = 0d0
*     Non-singlet
                  elseif(k.eq.3)then
*     C2
                     C2L(1) = C2ns1L
                     integC2(1) = C2ns1RS
*     CL
                     CLL(1) = CLns1L
                     integCL(1) = CLns1RS
*     C3
                     C3L(1) = C3ns1L
                     integC3(1) = C3ns1RS
                  endif
               endif
c*
c*     NNLO (Unknown yet)
c*
c               if(iptmx.ge.2)then
c*     C2
c                  C2L(2) = 0d0
c                  integC2(2) = 0d0
c*     CL
c                  CLL(2) = 0d0
c                  integCL(2) = 0d0
c*     C3
c                  C3L(2) = 0d0
c                  integC3(2) = 0d0
c               endif
*
*     Integrals
*
               do wipt=0,iptmx
                  SC2mCC(igrid,ixi,k,wipt,beta,alpha) = integC2(wipt) 
     1                 + C2L(wipt) * fL_CCm
                  SCLmCC(igrid,ixi,k,wipt,beta,alpha) = integCL(wipt)
     1                 + CLL(wipt) * fL_CCm 
                  SC3mCC(igrid,ixi,k,wipt,beta,alpha) = ( integC3(wipt)
     1                 + C3L(wipt) * fL_CCm ) * lambda
               enddo
            enddo
*
*     If the MSbar masses are to be used, add the appropriate term to
*     the NLO quark non-singlet coefficient functions.
*
            if(mass_scheme.eq."MSbar")then
               h1 = CF * 4d0
               if(MassRunning) h1 = h1 + CF * 3d0
     1              * dlog(xigrid(ixi)*krenQ)
               SC2mCC(igrid,ixi,3,1,beta,alpha) = 
     1              SC2mCC(igrid,ixi,3,1,beta,alpha)
     2              + 2d0 * h1 * ( 1d0 - lambda ) * dfL_CCm
               SCLmCC(igrid,ixi,3,1,beta,alpha) = 
     1              SCLmCC(igrid,ixi,3,1,beta,alpha)
     2              + 2d0 * h1 * ( 1d0 - lambda )**2 * dfL_CCm
               SC3mCC(igrid,ixi,3,1,beta,alpha) = 
     1              SC3mCC(igrid,ixi,3,1,beta,alpha)
     2              + 2d0 * h1 * ( 1d0 - lambda ) * ( dfL_CCm - fL_CCm )
     3              * lambda
            endif
         enddo
      endif
*
*     Massive zero FFNS needed for the FONLL scheme
*     (and also for the FFNS)
*
      if(MassScheme(1:4).eq."FFN0".or.MassScheme(1:5).eq."FONLL")then
*
*     Neutral Current
*
*     In case of FONLL-B
*
         ipt_FF = ipt
         if(MassScheme.eq."FONLL-B".and.ipt.ge.1) ipt_FF = 2
*
*     Variables needed for wrapping the integrand functions
*
         wnf    = Nf_FF
         walpha = alpha
         wbeta  = beta
*
*     Precompute integrals
*
         if(ipt_FF.ge.1)then
            wipt = 1
*
            sf = 1
            k  = 1
            wl = 1
            C2g1R0 = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            wl = 2
            C2g1RQ = dgauss(integrandsDISNCm0,a,b,eps(wipt))
*
            sf = 2
            k  = 1
            wl = 1
            CLg1R0 = dgauss(integrandsDISNCm0,a,b,eps(wipt))
*
            sf = 3
            k  = 3
            C3ns1RS = dgauss(integrandsDISzm,a,b,eps(wipt))
            C3ns1L  = C3NS1C(a)
         endif
         if(ipt_FF.ge.2)then
            wipt = 2
*
            sf = 1
            k  = 1
            wl = 1
            C2g2R0  = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            wl = 2
            C2g2RQ  = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            wl = 3
            C2g2RQ2 = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            wl = 4
            C2g2RF  = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            wl = 5
            C2g2RQF = dgauss(integrandsDISNCm0,a,b,eps(wipt))
*
            k  = 2
            wl = 1
            C2ps2R0  = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            wl = 2
            C2ps2RQ  = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            wl = 3
            C2ps2RQ2 = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            wl = 4
            C2ps2RF  = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            wl = 5
            C2ps2RQF = dgauss(integrandsDISNCm0,a,b,eps(wipt))
*
            k  = 3
            wl = 1
            C2ns2RS0  = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            C2ns2L0   = C2NS2CM0_A0(a)
            wl = 2
            C2ns2RSQ  = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            C2ns2LQ   = C2NS2CM0_AQ(a)
            wl = 3
            C2ns2RSQ2 = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            C2ns2LQ2  = C2NS2CM0_AQ2(a)
*
            sf = 2
            k  = 1
            wl = 1
            CLg2R0 = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            wl = 2
            CLg2RQ = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            wl = 4
            CLg2RF = dgauss(integrandsDISNCm0,a,b,eps(wipt))
*
            k = 2
            wl = 1
            CLps2R0 = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            wl = 2
            CLps2RQ = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            wl = 4
            CLps2RF = dgauss(integrandsDISNCm0,a,b,eps(wipt))
*
            k  = 3
            wl = 1
            CLns2R0 = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            wl = 2
            CLns2RQ = dgauss(integrandsDISNCm0,a,b,eps(wipt))
*
            sf = 3
            k  = 2
            C3ps2R  = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            k = 3
            C3ns2RS = dgauss(integrandsDISNCm0,a,b,eps(wipt))
            C3ns2L  = C3NSP2C(a,inf)
         endif
*
         do ixi=1,nxir
            lnQ  = dlog(xigrid(ixi*xistep))
            lnF  = lnkQF2 - lnQ
            lnQ2 = lnQ * lnQ
            lnQF = lnQ * lnF
            h1 = CF * 4d0
            if(MassRunning) h1 = h1
     1           + CF * 3d0 * dlog(xigrid(ixi*xistep)*krenQ)
            do k=1,3
               do wipt=1,ipt_FF
                  SC2m0NC(igrid,ixi,k,wipt,beta,alpha) = 0d0
                  SCLm0NC(igrid,ixi,k,wipt,beta,alpha) = 0d0
               enddo
            enddo
*
            do k=1,3
*
*     NLO
*
               if(ipt_FF.ge.1)then
*     Gluon
                  if(k.eq.1)then
*     C2
                     C2L(1) = 0d0
                     integC2(1) = C2g1R0 + C2g1RQ * lnQ
*     CL
                     integCL(1) = CLg1R0
*     Pure-Singlet
                  elseif(k.eq.2)then
*     C2
                     C2L(1) = 0d0
                     integC2(1) = 0d0
*     CL
                     integCL(1) = 0d0
*     Non-singlet
                  elseif(k.eq.3)then
*     C2
                     C2L(1) = 0d0
                     integC2(1) = 0d0
*     CL
                     integCL(1) = 0d0
                  endif
               endif
*
*     NNLO
*
               if(ipt_FF.ge.2)then
*     Gluon
                  if(k.eq.1)then
*     C2
                     C2L(2) = 0d0
                     integC2(2) = C2g2R0 + C2g2RQ * lnQ + C2g2RQ2 * lnQ2 
     1                          + C2g2RF * lnF + C2g2RQF * lnQF
*
*     If the MSbar masses are enabled, add the appropriate term
*
                     if(mass_scheme.eq."MSbar") 
     1                    integC2(2) = integC2(2) - 2d0 * h1 * C2g1RQ
*     CL
                     integCL(2) = CLg2R0 + CLg2RQ * lnQ + CLg2RF * lnF
*
*     Subtract the constant term for FONLL-B
*
                     if(MassScheme.eq."FONLL-B")then
                        integC2(2) = integC2(2) - C2g2R0
                        integCL(2) = integCL(2) - CLg2R0
                     endif
*     Pure-Singlet
                  elseif(k.eq.2)then
*     C2
                     C2L(2) = 0d0
                     integC2(2) = C2ps2R0 + C2ps2RQ * lnQ 
     1                          + C2ps2RQ2 * lnQ2 + C2ps2RF * lnF
     2                          + C2ps2RQF * lnQF
*     CL
                     integCL(2) = CLps2R0 + CLps2RQ * lnQ 
     1                          + CLps2RF * lnF
*
*     Subtract the constant term for FONLL-B
*
                     if(MassScheme.eq."FONLL-B")then
                        integC2(2) = integC2(2) - C2ps2R0
                        integCL(2) = integCL(2) - CLps2R0
                     endif
*     Non-singlet
                  elseif(k.eq.3)then
*     C2
                     C2L(2) = C2ns2L0 + C2ns2LQ * lnQ 
     1                      + C2ns2LQ2 * lnQ2
                     integC2(2) = C2ns2RS0 + C2ns2RSQ * lnQ 
     1                          + C2ns2RSQ2 * lnQ2
*     CL
                     integCL(2) = CLns2R0 + CLns2RQ * lnQ
*
*     Subtract the constant term for FONLL-B
*
                     if(MassScheme.eq."FONLL-B")then
                        integC2(2) = integC2(2) - C2ns2RS0
                        integCL(2) = integCL(2) - CLns2R0
                     endif
                  endif
               endif
*
*     Integrals
*
               do wipt=1,ipt_FF
                  SC2m0NC(igrid,ixi,k,wipt,beta,alpha) = integC2(wipt) 
     1               + C2L(wipt) * fL
                  SCLm0NC(igrid,ixi,k,wipt,beta,alpha) = integCL(wipt) 
               enddo
            enddo
         enddo
*
*     Charged Current
*
         walpha = alpha
         wbeta  = beta
*
         CGCOLM0 = dgauss(integrandsDISCCm0,a,b,eps(1))
*
         iptmx = ipt
         if(ipt.gt.1) iptmx = 1
*
         do ixi=1,nxir
            lnQ  = dlog(xigrid(ixi*xistep))
            lnC2 = lnkQF2 + lnQ
            lnC3 = lnkQF2 - lnQ
*
            do k=1,3
               do wipt=0,iptmx
                  SC2m0CC(igrid,ixi,k,wipt,beta,alpha) = 0d0
c                  SCLm0CC(igrid,ixi,k,wipt,beta,alpha) = 0d0
                  SC3m0CC(igrid,ixi,k,wipt,beta,alpha) = 0d0
               enddo
            enddo
*
            do wipt=0,iptmx
               CGM0 = 0d0
               if(wipt.ge.1) CGM0 = CGCOLM0
*
*     Gluon
*
               SC2m0CC(igrid,ixi,1,wipt,beta,alpha) = 
     1              SC2zm(igrid,Nf_FF,1,wipt,beta,alpha) + CGM0 * lnC2
               SCLm0CC(igrid,ixi,1,wipt,beta,alpha) = 
     1              SCLzm(igrid,Nf_FF,1,wipt,beta,alpha)
               SC3m0CC(igrid,ixi,1,wipt,beta,alpha) = 
     1              SC3zm(igrid,Nf_FF,1,wipt,beta,alpha) + CGM0 * lnC3
*
*     Non-singlet
*
               SC2m0CC(igrid,ixi,3,wipt,beta,alpha) = 
     1              SC2zm(igrid,Nf_FF,3,wipt,beta,alpha)
               SCLm0CC(igrid,ixi,3,wipt,beta,alpha) = 
     1              SCLzm(igrid,Nf_FF,3,wipt,beta,alpha)
               SC3m0CC(igrid,ixi,3,wipt,beta,alpha) = 
     1              SC3zm(igrid,Nf_FF,3,wipt,beta,alpha)
            enddo
         enddo
      endif
*
      return
      end
*
************************************************************************
*
*     Integrals needed for the computation of the target mass
*     corrections.
*
************************************************************************
      subroutine RSLintegralsDISTMC(beta,alpha)
*
      implicit none
*
      include "../commons/wrapDIS.h"
      include "../commons/grid.h"
      include "../commons/coeffhqmellin.h"
      include "../commons/integralsDIS.h"
**
*     Input Variables
*
      integer beta,alpha
**
*     Internal Variables
*
      integer bound
      double precision dgauss,c,d,eps
      double precision integrandsDISTMC
      external integrandsDISTMC
      data eps / 1d-5 /
*
      if(alpha.lt.beta) return
*
      bound = alpha-inter_degree(igrid)
      if(alpha.lt.inter_degree(igrid)) bound = 0
*
      c = max(xg(igrid,beta),xg(igrid,bound))
      d = min(1d0,xg(igrid,alpha+1))
*
      walpha = alpha
*
      J_TMC(igrid,beta,alpha) = dgauss(integrandsDISTMC,c,d,eps)
*
      return
      end
*
************************************************************************
*
*     Integrals of the small-x resummed coefficient functions.
*
************************************************************************
      subroutine RSLintegralsDISRes(beta,alpha,tau)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/gridAlpha.h"
      include "../commons/wrapResDIS.h"
      include "../commons/integralsResDIS.h"
      include "../commons/MassScheme.h"
**
*     Input Variables
*
      integer beta,alpha,tau
**
*     Internal Variables
*
      integer bound
      double precision dgauss,a,b,eps
      double precision integrandsDISzmRes
      external integrandsDISzmRes
      double precision integrandsDISNCcharmRes
      external integrandsDISNCcharmRes
      double precision integrandsDISNCcharm0Res
      external integrandsDISNCcharm0Res
      parameter(eps=1d-7)
*
*     Adjustment od the bounds of the integrals
*
      if(alpha.lt.beta) return
*
      bound = alpha-inter_degree(igrid)
      if(alpha.lt.inter_degree(igrid)) bound = 0
      a = max(xg(igrid,beta),xg(igrid,beta)/xg(igrid,alpha+1))
      b = min(1d0,xg(igrid,beta)/xg(igrid,bound))
*
*     Variables needed for wrapping the integrand functions
*
      walpha = alpha
      wbeta  = beta
      wtau   = tau
*
*     Precompute integrals
*
*     ZM-VFNS
*
      if(MassScheme.eq."ZM-VFNS".or.MassScheme(1:5).eq."FONLL".or.
     1     MassScheme(1:4).eq."FFNS".or.MassScheme(1:4).eq."FFN0")then
*
*     Initialize Integrals
*
         do k=1,2
            SC2zmRes(igrid,k,beta,alpha,tau) = 0d0
            SCLzmRes(igrid,k,beta,alpha,tau) = 0d0
         enddo
*     F2
         sf = 1
         do k=1,2
            SC2zmRes(igrid,k,beta,alpha,tau) =
     1           dgauss(integrandsDISzmRes,a,b,eps)
         enddo
*     FL
         sf = 2
         do k=1,2
            SCLzmRes(igrid,k,beta,alpha,tau) =
     1           dgauss(integrandsDISzmRes,a,b,eps)
         enddo
      endif
*
*     FFNS
*
      if(MassScheme(1:4).eq."FFNS".or.MassScheme(1:5).eq."FONLL")then
*
*     Initialize Integrals
*
         do k=1,2
            SC2charmNCRes(igrid,k,beta,alpha,tau) = 0d0
            SCLcharmNCRes(igrid,k,beta,alpha,tau) = 0d0
         enddo
*     F2
         sf = 1
         do k=1,2
            SC2charmNCRes(igrid,k,beta,alpha,tau) =
     1           dgauss(integrandsDISNCcharmRes,a,b,eps)
         enddo
*     FL
         sf = 2
         do k=1,2
            SCLcharmNCRes(igrid,k,beta,alpha,tau) =
     1           dgauss(integrandsDISNCcharmRes,a,b,eps)
         enddo
      endif
*
*     FFN0
*
      if(MassScheme(1:4).eq."FFN0".or.MassScheme(1:5).eq."FONLL")then
*
*     Initialize Integrals
*
         do k=1,2
            SC2charm0NCRes(igrid,k,beta,alpha,tau) = 0d0
            SCLcharm0NCRes(igrid,k,beta,alpha,tau) = 0d0
         enddo
*     F2
         sf = 1
         do k=1,2
            SC2charm0NCRes(igrid,k,beta,alpha,tau) =
     1           SC2zmRes(igrid,k,beta,alpha,tau)
     2           + dgauss(integrandsDISNCcharm0Res,a,b,eps)
         enddo
*     FL
         sf = 2
         do k=1,2
            SCLcharm0NCRes(igrid,k,beta,alpha,tau) =
     1           SCLzmRes(igrid,k,beta,alpha,tau)
c     2           + dgauss(integrandsDISNCcharm0Res,a,b,eps)
         enddo
      endif
*
      return
      end
