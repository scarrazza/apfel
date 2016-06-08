************************************************************************
*
*     RSLintegralsQCD.f:
*
*     This routine evaluates and dump on the common integrals.h once and 
*     for all for a given number of active flavours nf and for the pair 
*     (alpha,beta) the integral of the sum of Regular and Singular part 
*     plus the local part of all the splitting functions for all the needed
*     orders in QCD.
*     The index kk runs like that:
*
*     kk  =  1   2   3   4   5   6   7
*            +   -   V   qq  qg  gq  gg
* 
************************************************************************
      subroutine RSLintegralsQCD(nf,beta,alpha)
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/wrap.h"
      include "../commons/grid.h"
      include "../commons/integrals.h"
      include "../commons/kren.h"
**
*     Input Variables
*
      integer nf,beta,alpha
**
*     Internal Variables
*
      integer bound
      double precision PL(0:2),fL
      double precision X0NSC,X1NSC,P2NSPC,P2NSMC,X0GGC,X1GGC,P2GGC
      double precision dgauss,a,b,eps(0:2)
      double precision integrandsQCD
      double precision integ(0:2)
      double precision ns0L,gg0L,ns0RS,qg0R,gq0R,gg0RS
      double precision ns1L,gg1L,ns1RSp,ns1RSm,qq1RS,qg1RS,gq1RS,gg1RS
      double precision ns2Lp,ns2Lm,gg2L,ns2RSp,ns2RSm,ns2RSv
      double precision qq2RS,qg2RS,gq2RS,gg2RS
      double precision lnk,beta0apf,beta1apf
      external integrandsQCD
      data eps / 1d-6, 1d-5, 1d-4 /
*
*     Initialize Integrals
*
      do k=1,7
         do wipt=0,ipt
            SP(igrid,nf,k,wipt,beta,alpha) = 0d0
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
*     Variables needed for wrapping the integrand functions
*
      wnf    = nf
      walpha = alpha
      wbeta  = beta
*
*     Precompute integrals
*
      wipt = 0
      ns0L = X0NSC(a)
      gg0L = X0GGC(a,nf)
      k = 1
      ns0RS = dgauss(integrandsQCD,a,b,eps(wipt))
      k = 5
      qg0R  = dgauss(integrandsQCD,a,b,eps(wipt))
      k = 6
      gq0R  = dgauss(integrandsQCD,a,b,eps(wipt))
      k = 7
      gg0RS = dgauss(integrandsQCD,a,b,eps(wipt))
      if(ipt.ge.1)then
         wipt = 1
         ns1L = X1NSC(a,nf)
         gg1L = X1GGC(a,nf)
         k = 1
         ns1RSp = dgauss(integrandsQCD,a,b,eps(wipt))
         k = 2
         ns1RSm = dgauss(integrandsQCD,a,b,eps(wipt))
         k = 4
         qq1RS  = dgauss(integrandsQCD,a,b,eps(wipt))
         k = 5
         qg1RS  = dgauss(integrandsQCD,a,b,eps(wipt))
         k = 6
         gq1RS  = dgauss(integrandsQCD,a,b,eps(wipt))
         k = 7
         gg1RS  = dgauss(integrandsQCD,a,b,eps(wipt))
      endif
      if(ipt.ge.2)then
         wipt = 2
         ns2Lp = P2NSPC(a,nf)
         ns2Lm = P2NSMC(a,nf)
         gg2L  = P2GGC(a,nf)
         k = 1
         ns2RSp = dgauss(integrandsQCD,a,b,eps(wipt))
         k = 2
         ns2RSm = dgauss(integrandsQCD,a,b,eps(wipt))
         k = 3
         ns2RSv = dgauss(integrandsQCD,a,b,eps(wipt))
         k = 4
         qq2RS  = dgauss(integrandsQCD,a,b,eps(wipt))
         k = 5
         qg2RS  = dgauss(integrandsQCD,a,b,eps(wipt))
         k = 6
         gq2RS  = dgauss(integrandsQCD,a,b,eps(wipt))
         k = 7
         gg2RS  = dgauss(integrandsQCD,a,b,eps(wipt))
      endif
*
      do k=1,7
*
*     LO
*
*     Plus, Minus, Valence, Quark-Quark
         if(k.eq.1.or.k.eq.2.or.k.eq.3.or.k.eq.4)then
            PL(0)    = ns0L
            integ(0) = ns0RS
*     Quark-Gluon, Gluon-Quark
         elseif(k.eq.5)then
            PL(0)    = 0d0
            integ(0) = qg0R
         elseif(k.eq.6)then
            PL(0)    = 0d0
            integ(0) = gq0R
*     Gluon-Gluon
         elseif(k.eq.7)then
            PL(0)    = gg0L
            integ(0) = gg0RS
         endif
*
*     NLO
*
         if(ipt.ge.1)then
*     Plus
            if(k.eq.1)then
               PL(1)    = ns1L
               integ(1) = ns1RSp
*     Minus, Valence
            elseif(k.eq.2.or.k.eq.3)then
               PL(1)    = ns1L
               integ(1) = ns1RSm
*     Quark-Quark
            elseif(k.eq.4)then
               PL(1)    = ns1L
               integ(1) = qq1RS
*     Quark-Gluon
           elseif(k.eq.5)then
               PL(1)    = 0d0
               integ(1) = qg1RS
*     Gluon-Quark
           elseif(k.eq.6)then
               PL(1)    = 0d0
               integ(1) = gq1RS
*     Gluon-Gluon
            elseif(k.eq.7)then
               PL(1)    = gg1L
               integ(1) = gg1RS
            endif
         endif
*
*     NNLO
*
         if(ipt.ge.2)then
*     Plus
            if(k.eq.1)then
               PL(2)    = ns2Lp
               integ(2) = ns2RSp
*     Minus
            elseif(k.eq.2)then
               PL(2)    = ns2Lm
               integ(2) = ns2RSm
*     Valence
            elseif(k.eq.3)then
               PL(2)    = ns2Lm
               integ(2) = ns2RSv
*     Quark-Quark
            elseif(k.eq.4)then
               PL(2)    = ns2Lp
               integ(2) = qq2RS
*     Quark-Gluon
            elseif(k.eq.5)then
               PL(2)    = 0d0
               integ(2) = qg2RS
*     Gluon-Quark
            elseif(k.eq.6)then
               PL(2)    = 0d0
               integ(2) = gq2RS
*     Gluon-Gluon
            elseif(k.eq.7)then
               PL(2)    = gg2L
               integ(2) = gg2RS
            endif
         endif
*
*     Integrals
*
         do wipt=0,ipt
            SP(igrid,nf,k,wipt,beta,alpha) = integ(wipt) + PL(wipt) * fL
         enddo
*
*     In case of muR.ne.muF...
*     Eq. (2.8) of hep-ph/0408244
*
         if(kren.ne.1d0)then
            lnk = - dlog(kren)
            if(ipt.eq.1)then
               SP(igrid,nf,k,1,beta,alpha) = SP(igrid,nf,k,1,beta,alpha) 
     1         - beta0apf(nf) * lnk * SP(igrid,nf,k,0,beta,alpha) 
            elseif(ipt.eq.2)then
               SP(igrid,nf,k,2,beta,alpha) = SP(igrid,nf,k,2,beta,alpha) 
     1         - 2d0 * beta0apf(nf) * lnk * SP(igrid,nf,k,1,beta,alpha)
     2         - ( beta1apf(nf) - beta0apf(nf)**2 * lnk )
     3         * lnk * SP(igrid,nf,k,0,beta,alpha)
               SP(igrid,nf,k,1,beta,alpha) = SP(igrid,nf,k,1,beta,alpha) 
     1         - beta0apf(nf) * lnk * SP(igrid,nf,k,0,beta,alpha) 
            endif
         endif
      enddo
*
      return
      end
*
************************************************************************
*
*     Integrals for the time-like evolution.
*
************************************************************************
      subroutine RSLintegralsQCDT(nf,beta,alpha)
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/wrap.h"
      include "../commons/grid.h"
      include "../commons/integrals.h"
      include "../commons/kren.h"
**
*     Input Variables
*
      integer nf,beta,alpha
**
*     Internal Variables
*
      integer bound
      double precision PL(0:2),fL
      double precision X0NSC,X1NSTC,P2NSPTC,P2NSMTC,X0GGC,X1GGTC,P2GGTC
      double precision dgauss,a,b,eps(0:2)
      double precision integrandsQCDT
      double precision integ(0:2)
      double precision ns0L,gg0L,ns0RS,qg0R,gq0R,gg0RS
      double precision ns1L,gg1L,ns1RSp,ns1RSm,qq1RS,qg1RS,gq1RS,gg1RS
      double precision ns2Lp,ns2Lm,gg2L,ns2RSp,ns2RSm,ns2RSv
      double precision qq2RS,qg2RS,gq2RS,gg2RS
      double precision lnk,beta0apf,beta1apf
      external integrandsQCDT
      data eps / 1d-7, 1d-6, 1d-5 /
*
*     Initialize Integrals
*
      do k=1,7
         do wipt=0,ipt
            SP(igrid,nf,k,wipt,beta,alpha) = 0d0
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
*     Variables needed for wrapping the integrand functions
*
      wnf    = nf
      walpha = alpha
      wbeta  = beta
*
*     Precompute integrals
*
      wipt = 0
      ns0L = X0NSC(a)
      gg0L = X0GGC(a,nf)
      k = 1
      ns0RS = dgauss(integrandsQCDT,a,b,eps(wipt))
      k = 5
      qg0R  = dgauss(integrandsQCDT,a,b,eps(wipt))
      k = 6
      gq0R  = dgauss(integrandsQCDT,a,b,eps(wipt))
      k = 7
      gg0RS = dgauss(integrandsQCDT,a,b,eps(wipt))
      if(ipt.ge.1)then
         wipt = 1
         ns1L = X1NSTC(a,nf)
         gg1L = X1GGTC(a,nf)
         k = 1
         ns1RSp = dgauss(integrandsQCDT,a,b,eps(wipt))
         k = 2
         ns1RSm = dgauss(integrandsQCDT,a,b,eps(wipt))
         k = 4
         qq1RS  = dgauss(integrandsQCDT,a,b,eps(wipt))
         k = 5
         qg1RS  = dgauss(integrandsQCDT,a,b,eps(wipt))
         k = 6
         gq1RS  = dgauss(integrandsQCDT,a,b,eps(wipt))
         k = 7
         gg1RS  = dgauss(integrandsQCDT,a,b,eps(wipt))
      endif
      if(ipt.ge.2)then
         wipt = 2
         ns2Lp = P2NSPTC(a,nf)
         ns2Lm = P2NSMTC(a,nf)
         gg2L  = P2GGTC(a,nf)
         k = 1
         ns2RSp = dgauss(integrandsQCDT,a,b,eps(wipt))
         k = 2
         ns2RSm = dgauss(integrandsQCDT,a,b,eps(wipt))
         k = 3
         ns2RSv = dgauss(integrandsQCDT,a,b,eps(wipt))
         k = 4
         qq2RS  = dgauss(integrandsQCDT,a,b,eps(wipt))
         k = 5
         qg2RS  = dgauss(integrandsQCDT,a,b,eps(wipt))
         k = 6
         gq2RS  = dgauss(integrandsQCDT,a,b,eps(wipt))
         k = 7
         gg2RS  = dgauss(integrandsQCDT,a,b,eps(wipt))
      endif
*
      do k=1,7
*
*     LO
*
*     Plus, Minus, Valence, Quark-Quark
         if(k.eq.1.or.k.eq.2.or.k.eq.3.or.k.eq.4)then
            PL(0)    = ns0L
            integ(0) = ns0RS
*     Quark-Gluon, Gluon-Quark
         elseif(k.eq.6)then
            PL(0)    = 0d0
            integ(0) = qg0R / 2d0 / nf
         elseif(k.eq.5)then
            PL(0)    = 0d0
            integ(0) = 2d0 * nf * gq0R
*     Gluon-Gluon
         elseif(k.eq.7)then
            PL(0)    = gg0L
            integ(0) = gg0RS
         endif
*
*     NLO
*
         if(ipt.ge.1)then
*     Plus
            if(k.eq.1)then
               PL(1)    = ns1L
               integ(1) = ns1RSp
*     Minus, Valence
            elseif(k.eq.2.or.k.eq.3)then
               PL(1)    = ns1L
               integ(1) = ns1RSm
*     Quark-Quark
            elseif(k.eq.4)then
               PL(1)    = ns1L
               integ(1) = qq1RS
*     Quark-Gluon
           elseif(k.eq.6)then
               PL(1)    = 0d0
               integ(1) = qg1RS
*     Gluon-Quark
           elseif(k.eq.5)then
               PL(1)    = 0d0
               integ(1) = gq1RS
*     Gluon-Gluon
            elseif(k.eq.7)then
               PL(1)    = gg1L
               integ(1) = gg1RS
            endif
         endif
*
*     NNLO
*
         if(ipt.ge.2)then
*     Plus
            if(k.eq.1)then
               PL(2)    = ns2Lp
               integ(2) = ns2RSp
*     Minus
            elseif(k.eq.2)then
               PL(2)    = ns2Lm
               integ(2) = ns2RSm
*     Valence
            elseif(k.eq.3)then
               PL(2)    = ns2Lm
               integ(2) = ns2RSv
*     Quark-Quark
            elseif(k.eq.4)then
               PL(2)    = ns2Lp
               integ(2) = qq2RS
*     Quark-Gluon
            elseif(k.eq.6)then
               PL(2)    = 0d0
               integ(2) = qg2RS
*     Gluon-Quark
            elseif(k.eq.5)then
               PL(2)    = 0d0
               integ(2) = gq2RS
*     Gluon-Gluon
            elseif(k.eq.7)then
               PL(2)    = gg2L
               integ(2) = gg2RS
            endif
         endif
*
*     Integrals
*
         do wipt=0,ipt
            SP(igrid,nf,k,wipt,beta,alpha) = integ(wipt) + PL(wipt) * fL
         enddo
*
*     In case of muR.ne.muF...
*     Eq. (2.8) of hep-ph/0408244
*
         if(kren.ne.1d0)then
            lnk = - dlog(kren)
            if(ipt.eq.1)then
               SP(igrid,nf,k,1,beta,alpha) = SP(igrid,nf,k,1,beta,alpha) 
     1         - beta0apf(nf) * lnk * SP(igrid,nf,k,0,beta,alpha) 
            elseif(ipt.eq.2)then
               SP(igrid,nf,k,2,beta,alpha) = SP(igrid,nf,k,2,beta,alpha) 
     1         - 2d0 * beta0apf(nf) * lnk * SP(igrid,nf,k,1,beta,alpha)
     2         - ( beta1apf(nf) - beta0apf(nf)**2 * lnk )
     3         * lnk * SP(igrid,nf,k,0,beta,alpha)
               SP(igrid,nf,k,1,beta,alpha) = SP(igrid,nf,k,1,beta,alpha) 
     1         - beta0apf(nf) * lnk * SP(igrid,nf,k,0,beta,alpha) 
            endif
         endif
      enddo
*
      return
      end
*
************************************************************************
*
*     Integrals of the small-x resummed splitting functions.
*
************************************************************************
      subroutine RSLintegralsQCDRes(la,beta,alpha,tau)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/gridAlpha.h"
      include "../commons/wrapRes.h"
      include "../commons/integralsRes.h"
**
*     Input Variables
*
      integer la,beta,alpha,tau
**
*     Internal Variables
*
      integer bound
      double precision dgauss,a,b,eps
      double precision integrandsQCDRes
      external integrandsQCDRes
      parameter(eps=1d-7)
*
*     Initialize Integrals
*
      do k=4,7
         SPRes(igrid,k,la,beta,alpha,tau) = 0d0
      enddo
*
*     Adjustment of the bounds of the integrals
*
      if(alpha.lt.beta)then
         return
      else
         bound = alpha - inter_degree(igrid)
         if(alpha.lt.inter_degree(igrid)) bound = 0
         a = max(xg(igrid,beta),xg(igrid,beta)/xg(igrid,alpha+1))
         b = min(1d0,xg(igrid,beta)/xg(igrid,bound))
      endif
*
*     Variables needed for wrapping the integrand functions
*
      walpha = alpha
      wbeta  = beta
      wtau   = tau
*
*     Precompute integrals
*
      do k=4,7
         SPRes(igrid,k,la,beta,alpha,tau) =
     1        dgauss(integrandsQCDRes,a,b,eps)
      enddo
*
      return
      end
*
************************************************************************
*
*     Integrals for the polarized evolution.
*
************************************************************************
      subroutine RSLintegralsQCDPol(nf,beta,alpha)
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/wrap.h"
      include "../commons/grid.h"
      include "../commons/integrals.h"
      include "../commons/kren.h"
**
*     Input Variables
*
      integer nf,beta,alpha
**
*     Internal Variables
*
      integer bound
      double precision PL(0:2),fL
      double precision X0NSPC,X1NSPC,P2NSPPC,P2NSMPC,X0GGPC,X1GGPC
      double precision P2GGPC
      double precision dgauss,a,b,eps(0:2)
      double precision integrandsQCDPol
      double precision integ(0:2)
      double precision ns0L,gg0L,ns0RS,qg0R,gq0R,gg0RS
      double precision ns1L,gg1L,ns1RSp,ns1RSm,qq1RS,qg1RS,gq1RS,gg1RS
      double precision ns2Lp,ns2Lm,gg2L,ns2RSp,ns2RSm,ns2RSv
      double precision qq2RS,qg2RS,gq2RS,gg2RS
      double precision lnk,beta0apf,beta1apf
      external integrandsQCDPol
      data eps / 1d-6, 1d-5, 1d-4 /
*
*     Initialize Integrals
*
      do k=1,7
         do wipt=0,ipt
            SP(igrid,nf,k,wipt,beta,alpha) = 0d0
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
*     Variables needed for wrapping the integrand functions
*
      wnf    = nf
      walpha = alpha
      wbeta  = beta
*
*     Precompute integrals
*
      wipt = 0
      ns0L = X0NSPC(a)
      gg0L = X0GGPC(a,nf)
      k = 1
      ns0RS = dgauss(integrandsQCDPol,a,b,eps(wipt))
      k = 5
      qg0R  = dgauss(integrandsQCDPol,a,b,eps(wipt))
      k = 6
      gq0R  = dgauss(integrandsQCDPol,a,b,eps(wipt))
      k = 7
      gg0RS = dgauss(integrandsQCDPol,a,b,eps(wipt))
      if(ipt.ge.1)then
         wipt = 1
         ns1L = X1NSPC(a,nf)
         gg1L = X1GGPC(a,nf)
         k = 1
         ns1RSp = dgauss(integrandsQCDPol,a,b,eps(wipt))
         k = 2
         ns1RSm = dgauss(integrandsQCDPol,a,b,eps(wipt))
         k = 4
         qq1RS  = dgauss(integrandsQCDPol,a,b,eps(wipt))
         k = 5
         qg1RS  = dgauss(integrandsQCDPol,a,b,eps(wipt))
         k = 6
         gq1RS  = dgauss(integrandsQCDPol,a,b,eps(wipt))
         k = 7
         gg1RS  = dgauss(integrandsQCDPol,a,b,eps(wipt))
      endif
      if(ipt.ge.2)then
         wipt = 2
         ns2Lp = P2NSPPC(a,nf)
         ns2Lm = P2NSMPC(a,nf)
         gg2L  = P2GGPC(a,nf)
         k = 1
         ns2RSp = dgauss(integrandsQCDPol,a,b,eps(wipt))
         k = 2
         ns2RSm = dgauss(integrandsQCDPol,a,b,eps(wipt))
         k = 3
         ns2RSv = dgauss(integrandsQCDPol,a,b,eps(wipt))
         k = 4
         qq2RS  = dgauss(integrandsQCDPol,a,b,eps(wipt))
         k = 5
         qg2RS  = dgauss(integrandsQCDPol,a,b,eps(wipt))
         k = 6
         gq2RS  = dgauss(integrandsQCDPol,a,b,eps(wipt))
         k = 7
         gg2RS  = dgauss(integrandsQCDPol,a,b,eps(wipt))
      endif
*
      do k=1,7
*
*     LO
*
*     Plus, Minus, Valence, Quark-Quark
         if(k.eq.1.or.k.eq.2.or.k.eq.3.or.k.eq.4)then
            PL(0)    = ns0L
            integ(0) = ns0RS
*     Quark-Gluon, Gluon-Quark
         elseif(k.eq.5)then
            PL(0)    = 0d0
            integ(0) = qg0R
         elseif(k.eq.6)then
            PL(0)    = 0d0
            integ(0) = gq0R
*     Gluon-Gluon
         elseif(k.eq.7)then
            PL(0)    = gg0L
            integ(0) = gg0RS
         endif
*
*     NLO
*
         if(ipt.ge.1)then
*     Plus
            if(k.eq.1)then
               PL(1)    = ns1L
               integ(1) = ns1RSp
*     Minus, Valence
            elseif(k.eq.2.or.k.eq.3)then
               PL(1)    = ns1L
               integ(1) = ns1RSm
*     Quark-Quark
            elseif(k.eq.4)then
               PL(1)    = ns1L
               integ(1) = qq1RS
*     Quark-Gluon
           elseif(k.eq.5)then
               PL(1)    = 0d0
               integ(1) = qg1RS
*     Gluon-Quark
           elseif(k.eq.6)then
               PL(1)    = 0d0
               integ(1) = gq1RS
*     Gluon-Gluon
            elseif(k.eq.7)then
               PL(1)    = gg1L
               integ(1) = gg1RS
            endif
         endif
*
*     NNLO
*
         if(ipt.ge.2)then
*     Plus
            if(k.eq.1)then
               PL(2)  = ns2Lp
               integ(2) = ns2RSp
*     Minus
            elseif(k.eq.2)then
               PL(2)    = ns2Lm
               integ(2) = ns2RSm
*     Valence
            elseif(k.eq.3)then
               PL(2)    = ns2Lm
               integ(2) = ns2RSv
*     Quark-Quark
            elseif(k.eq.4)then
               PL(2)    = ns2Lp
               integ(2) = qq2RS
*     Quark-Gluon
            elseif(k.eq.5)then
               PL(2)    = 0d0
               integ(2) = qg2RS
*     Gluon-Quark
            elseif(k.eq.6)then
               PL(2)    = 0d0
               integ(2) = gq2RS
*     Gluon-Gluon
            elseif(k.eq.7)then
               PL(2)    = gg2L
               integ(2) = gg2RS
            endif
         endif
*
*     Integrals
*
         do wipt=0,ipt
            SP(igrid,nf,k,wipt,beta,alpha) = integ(wipt) + PL(wipt) * fL
         enddo
*
*     In case of muR.ne.muF...
*     Eq. (2.8) of hep-ph/0408244
*
         if(kren.ne.1d0)then
            lnk = - dlog(kren)
            if(ipt.eq.1)then
               SP(igrid,nf,k,1,beta,alpha) = SP(igrid,nf,k,1,beta,alpha) 
     1         - beta0apf(nf) * lnk * SP(igrid,nf,k,0,beta,alpha) 
            elseif(ipt.eq.2)then
               SP(igrid,nf,k,2,beta,alpha) = SP(igrid,nf,k,2,beta,alpha) 
     1         - 2d0 * beta0apf(nf) * lnk * SP(igrid,nf,k,1,beta,alpha)
     2         - ( beta1apf(nf) - beta0apf(nf)**2 * lnk )
     3         * lnk * SP(igrid,nf,k,0,beta,alpha)
               SP(igrid,nf,k,1,beta,alpha) = SP(igrid,nf,k,1,beta,alpha) 
     1         - beta0apf(nf) * lnk * SP(igrid,nf,k,0,beta,alpha) 
            endif
         endif
      enddo
*
      return
      end
