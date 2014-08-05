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
      subroutine RSLintegralsMatching(beta,alpha)
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
      integer beta,alpha
**
*     Internal Variables
*
      integer bound
      double precision ML(5,0:2),fL
      double precision ANS2qqH_L,AS2ggH_L
      double precision dgauss,a,b,eps
      double precision integrandsMatching
      external integrandsMatching
      parameter(eps=1d-4)
*
*     Initialize Integrals
*
      do k=1,5
         do wipt=0,ipt
            SM(igrid,k,wipt,beta,alpha) = 0d0
         enddo
      enddo
      if(alpha.eq.beta)then
         SM(igrid,1,0,beta,alpha) = 1d0
         SM(igrid,2,0,beta,alpha) = 1d0
         SM(igrid,5,0,beta,alpha) = 1d0
      endif
      if(ipt.le.1) return
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
      walpha = alpha
      wbeta  = beta
      wipt   = 2
*
      do k=1,5
*
*     Local contributions
*
*
*     NNLO
*
*     Valence, Quark-Quark
         if(k.eq.1.or.k.eq.2)then
            ML(k,2) = ANS2qqH_L(a)
*     Quark-Gluon, Gluon-Quark
         elseif(k.eq.3.or.k.eq.4)then
            ML(k,2) = 0d0
*     Gluon-Gluon
         elseif(k.eq.5)then
            ML(k,2) = AS2ggH_L(a)
         endif
*
*     Integrals
*
         SM(igrid,k,wipt,beta,alpha)= dgauss(integrandsMatching,a,b,eps) 
     1                              + ML(k,wipt) * fL
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
      subroutine RSLintegralsMatchingT(beta,alpha)
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
      integer beta,alpha
**
*     Internal Variables
*
      integer bound
      double precision ML(5,0:2)
      double precision ANS2qqH_L,AS2ggH_L
      double precision dgauss,a,b,eps
      double precision integrandsMatchingT
      external integrandsMatchingT
      parameter(eps=1d-4)
*
*     Initialize Integrals
*
      do k=1,5
         do wipt=0,ipt
            SM(igrid,k,wipt,beta,alpha) = 0d0
         enddo
      enddo
      if(alpha.eq.beta)then
         SM(igrid,1,0,beta,alpha) = 1d0
         SM(igrid,2,0,beta,alpha) = 1d0
         SM(igrid,5,0,beta,alpha) = 1d0
      endif
      if(ipt.lt.1) return
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
      walpha = alpha
      wbeta  = beta
      wipt   = 1
*
*     Integrals
*
c      do k=1,5
      do k=3,3
         SM(igrid,k,wipt,beta,alpha) =
     1        dgauss(integrandsMatchingT,a,b,eps) 
      enddo
*
      return
      end
