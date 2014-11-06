************************************************************************
*
*     integrandsMatching.f:
*
*     This function returns the integrands need to match PDF at the
*     heavy quark thresholds.
*     In order to wrap it in a form that can be fed to the integration 
*     function DGAUSS, it is an explicit function only of the integration 
*     variable y while it depends on the other variables:
*
*     1) the perturbatibe order pt, 
*     2) the grid indices walpha and wbeta,
*     3) the particular plitting function denoted by k such that:
*
*        k  =  1   2   3   4   5
*              V   qq  qg  gq  gg
*
*     that are contained in the common block wrap.h
* 
************************************************************************
      function integrandsMatching(y)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/wrap.h"
**
*     Input Variables
*
      double precision y
**
*     Internal Variables
*
      double precision z,w_int,fR,fS,fL
      double precision MR(7,0:2),MS(7,0:2)
      double precision ANS2qqH_R,ANS2qqH_S,APS2Hq,AS2Hg,AS2gqH
      double precision AS2ggH_R,AS2ggH_S
**
*     Output Variables
*
      double precision integrandsMatching
*
*     Interpolant functions
*
      z = xg(igrid,wbeta) / y
*
      fL = 0d0
      if(walpha.eq.wbeta) fL = 1d0
*
      fR = w_int(inter_degree(igrid),walpha,z)
      fS = fR - fL
*
*     Contructing integrands order by order
*
*     NLO (always zero for the space-like evolution)
*
      if(wipt.ge.1)then
         MR(k,1) = 0d0
         MS(k,1) = 0d0
      endif
*
*     NNLO
*
      if(wipt.ge.2)then
*     Valence
         if(k.eq.1)then
            MR(k,2) = ANS2qqH_R(y)
            MS(k,2) = ANS2qqH_S(y)
*     Quark-Quark
         elseif(k.eq.2)then
            MR(k,2) = ANS2qqH_R(y) + APS2Hq(y)
            MS(k,2) = ANS2qqH_S(y)
*     Quark-Gluon
         elseif(k.eq.3)then
            MR(k,2) = AS2Hg(y)
            MS(k,2) = 0d0
*     Gluon-Quark
         elseif(k.eq.4)then
            MR(k,2) = AS2gqH(y)
            MS(k,2) = 0d0
*     Gluon-Gluon
         elseif(k.eq.5)then
            MR(k,2) = AS2ggH_R(y)
            MS(k,2) = AS2ggH_S(y)
         endif
      endif
*
      integrandsMatching = MR(k,wipt) * fR + MS(k,wipt) * fS
c      integrandsMatching = MR(k,wipt) * fR
*
      return
      end
*
************************************************************************
*
*     Integrands for the time-like evolution.
*
************************************************************************      
      function integrandsMatchingT(y)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/wrap.h"
**
*     Input Variables
*
      double precision y
**
*     Internal Variables
*
      double precision z,w_int,fR!,fS,fL
      double precision MR(7,0:2)!,MS(7,0:2)
      double precision AS1HgT
**
*     Output Variables
*
      double precision integrandsMatchingT
*
*     Interpolant functions
*
      z = xg(igrid,wbeta) / y
c*
c      fL = 0d0
c      if(walpha.eq.wbeta) fL = 1d0
*
      fR = w_int(inter_degree(igrid),walpha,z)
c      fS = fR - fL
*
*     Contructing integrands order by order
*
*     NLO
*
      if(wipt.ge.1)then
*     Quark-Gluon
         if(k.eq.3)then
            MR(k,1) = AS1HgT(y)
         else
            MR(k,1) = 0d0
         endif
c         MS(k,1) = 0d0
      endif
*
c      integrandsMatchingT = MR(k,wipt) * fR + MS(k,wipt) * fS
      integrandsMatchingT = MR(k,wipt) * fR
*
      return
      end
