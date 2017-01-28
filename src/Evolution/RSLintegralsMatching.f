************************************************************************
*
*     RSLintegralsMatching.f:
*
*     This routine evaluates and dump on the common integrals.h once and 
*     for all for a given number of active flavours for the pair 
*     (alpha,beta) the integral of the sum of Regular and Singular part 
*     plus the local part of the matching conditions.
*     This routine should not be called for the LO order evolution but
*     in case it is it returns the "identity" matching.
*     The index kk runs like that:
*
*     kk  =  1   2   3   4   5 
*            V   qq  qg  gq  gg
* 
************************************************************************
      subroutine RSLintegralsMatching(nf,beta,alpha)
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/wrap.h"
      include "../commons/grid.h"
      include "../commons/integrals.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      integer nf,beta,alpha
**
*     Internal Variables
*
      integer bound
      double precision ML(2),fL
      double precision ANS2qqH_L,AS2ggH_L
      double precision AS1ggH_mass_L,ANS2qqH_mass_L,AS2ggH_mass_L
      double precision dgauss,a,b,eps
      double precision integrandsMatching
      external integrandsMatching
      parameter(eps=1d-4)
*
*     Initialize Integrals
*
      do k=1,5
         do wipt=0,ipt
            SM(igrid,nf,k,wipt,beta,alpha) = 0d0
         enddo
      enddo
      if(alpha.eq.beta)then
         SM(igrid,nf,1,0,beta,alpha) = 1d0
         SM(igrid,nf,2,0,beta,alpha) = 1d0
         SM(igrid,nf,5,0,beta,alpha) = 1d0
      endif
      if(ipt.eq.0) return
      if(ipt.eq.1.and.k2th(nf).eq.1d0) return
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
*     Local contributions
*
      do k=1,5
*
*     NLO
*
         wipt = 1
*     Gluon-Gluon
         if(k.eq.5.and.k2th(nf).ne.1d0)then
            ML(1) = AS1ggH_mass_L(wnf,a)
         else
            ML(1) = 0d0
         endif
*
*     Integrals
*
         if(k2th(nf).ne.1d0.and.(k.eq.3.or.k.eq.5))then
            SM(igrid,nf,k,1,beta,alpha) =
     1           dgauss(integrandsMatching,a,b,eps)
     2           + ML(1) * fL
         endif
*
*     NNLO
*
         if(ipt.ge.2)then
            wipt = 2
*     Valence, Quark-Quark
            if(k.eq.1.or.k.eq.2)then
               ML(2) = ANS2qqH_L(a)
*     Quark-Gluon, Gluon-Quark
            elseif(k.eq.3.or.k.eq.4)then
               ML(2) = 0d0
*     Gluon-Gluon
            elseif(k.eq.5)then
               ML(2) = AS2ggH_L(a)
            endif
            if(k2th(nf).ne.1d0)then
               if(k.eq.1.or.k.eq.2)then
                  ML(2) = ML(2) + ANS2qqH_mass_L(wnf,a)
*     Gluon-Gluon
               elseif(k.eq.5)then
                  ML(2) = ML(2) + AS2ggH_mass_L(wnf,a)
               endif
            endif
*
            SM(igrid,nf,k,2,beta,alpha) =
     1           dgauss(integrandsMatching,a,b,eps)
     2           + ML(2) * fL
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
      subroutine RSLintegralsMatchingT(nf,beta,alpha)
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/wrap.h"
      include "../commons/grid.h"
      include "../commons/integrals.h"
**
*     Input Variables
*
      integer nf,beta,alpha
**
*     Internal Variables
*
      integer bound
      double precision dgauss,a,b,eps
      double precision integrandsMatchingT
      external integrandsMatchingT
      parameter(eps=1d-5)
*
*     Initialize Integrals
*
      do k=1,5
         do wipt=0,ipt
            SM(igrid,nf,k,wipt,beta,alpha) = 0d0
         enddo
      enddo
      if(alpha.eq.beta)then
         SM(igrid,nf,1,0,beta,alpha) = 1d0
         SM(igrid,nf,2,0,beta,alpha) = 1d0
         SM(igrid,nf,5,0,beta,alpha) = 1d0
      endif
      if(ipt.eq.0) return
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
*     Variables needed for wrapping the integrand functions
*
      wnf    = nf
      walpha = alpha
      wbeta  = beta
      wipt   = 1
*
*     Integrals
*
      do k=3,3
         SM(igrid,nf,k,wipt,beta,alpha) =
     1        dgauss(integrandsMatchingT,a,b,eps) 
      enddo
*
      return
      end
*
************************************************************************
*
*     Integrals of the small-x resummed matching conditions.
*
************************************************************************
      subroutine RSLintegralsMatchingRes(nf,beta,alpha)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/gridAlpha.h"
      include "../commons/wrapRes.h"
      include "../commons/integrals.h"
**
*     Input Variables
*
      integer nf,beta,alpha
**
*     Internal Variables
*
      integer bound
      double precision dgauss,a,b,eps
      double precision integrandsMatchingRes
      external integrandsMatchingRes
      parameter(eps=1d-7)
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
      wnf    = nf
      walpha = alpha
      wbeta  = beta
*
*     Precompute integrals and put them in the LO slot of
*     the matching condition array because they don't have to be
*     multiplied by alpha_s.
*
*     Kqq
      k = 2
      SM(igrid,nf,2,0,beta,alpha) = SM(igrid,nf,2,0,beta,alpha)
     1     + dgauss(integrandsMatchingRes,a,b,eps)
*     Kqg
      k = 1
      SM(igrid,nf,3,0,beta,alpha) = SM(igrid,nf,3,0,beta,alpha)
     1     + dgauss(integrandsMatchingRes,a,b,eps)
*
      return
      end
