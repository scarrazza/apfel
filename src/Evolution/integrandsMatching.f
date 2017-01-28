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
      include "../commons/m2th.h"
**
*     Input Variables
*
      double precision y
**
*     Internal Variables
*
      double precision z,w_int,fR,fS,fL
      double precision MR,MS
      double precision ANS2qqH_R,ANS2qqH_S,APS2Hq,AS2Hg,AS2gqH
      double precision AS2ggH_R,AS2ggH_S
      double precision AS1Hg_mass
      double precision APS2Hq_mass,AS2Hg_mass
      double precision ANS2qqH_mass_R,ANS2qqH_mass_S
      double precision AS2gqH_mass
      double precision AS2ggH_mass_R,AS2ggH_mass_S
**
*     Output Variables
*
      double precision integrandsMatching
*
*     LO always zero
*
      if(wipt.eq.0)then
         integrandsMatching  = 0d0
         return
      endif
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
*     NLO
*
      if(wipt.eq.1)then
*     Log terms
         if(k2th(wnf).ne.1d0.and.k.eq.3)then
*     Quark-Gluon
            MR = AS1Hg_mass(wnf,y)
         else
            MR = 0d0
         endif
         MS = 0d0
*
*     NNLO
*
      elseif(wipt.ge.2)then
*     Valence
         if(k.eq.1)then
            MR = ANS2qqH_R(y)
            MS = ANS2qqH_S(y)
*     Quark-Quark
         elseif(k.eq.2)then
            MR = ANS2qqH_R(y) + APS2Hq(y)
            MS = ANS2qqH_S(y)
*     Quark-Gluon
         elseif(k.eq.3)then
            MR = AS2Hg(y)
            MS = 0d0
*     Gluon-Quark
         elseif(k.eq.4)then
            MR = AS2gqH(y)
            MS = 0d0
*     Gluon-Gluon
         elseif(k.eq.5)then
            MR = AS2ggH_R(y)
            MS = AS2ggH_S(y)
         endif
*     Log terms
         if(k2th(wnf).ne.1d0)then
*     Valence
            if(k.eq.1)then
               MR = MR + ANS2qqH_mass_R(wnf,y)
               MS = MS + ANS2qqH_mass_S(wnf,y)
*     Quark-Quark
            elseif(k.eq.2)then
               MR = MR + ANS2qqH_mass_R(wnf,y) + APS2Hq_mass(wnf,y)
               MS = MS + ANS2qqH_mass_S(wnf,y)
*     Quark-Gluon
            elseif(k.eq.3)then
               MR = MR + AS2Hg_mass(wnf,y)
*     Gluon-Quark
            elseif(k.eq.4)then
               MR = MR + AS2gqH_mass(wnf,y)
*     Gluon-Gluon
            elseif(k.eq.5)then
               MR = MR + AS2ggH_mass_R(wnf,y)
               MS = MS + AS2ggH_mass_S(wnf,y)
            endif
         endif
      endif
*
      integrandsMatching = MR * fR + MS * fS
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
      include "../commons/m2th.h"
**
*     Input Variables
*
      double precision y
**
*     Internal Variables
*
      double precision z,w_int,fR!,fS,fL
      double precision MR!,MS
      double precision AS1HgT,AS1HgT_mass
**
*     Output Variables
*
      double precision integrandsMatchingT
*
*     At LO always zero.
*     At NNLO they are not known and we set them to zero.
*
      if(wipt.eq.0.or.wipt.eq.2)then
         integrandsMatchingT = 0d0
         return
      endif
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
      if(wipt.eq.1)then
*     Quark-Gluon
         if(k.eq.3)then
            MR = AS1HgT(y)
*     Log terms
            if(k2th(wnf).ne.1d0)then
               MR = MR + AS1HgT_mass(wnf,y)
            endif
         else
            MR = 0d0
         endif
c         MS = 0d0
      endif
*
c      integrandsMatchingT = MR(k,wipt) * fR + MS(k,wipt) * fS
      integrandsMatchingT = MR * fR
*
      return
      end
*
************************************************************************
*
*     Integrands for the small-x resummed matching conditions.
*     The following routine calls the function "xDeltaK" that is the
*     Fortran wrapper of the c++ function provided by the Bonvini's
*     code HELL.
*     The input variables of the function "xDeltaP" are:
*
*        xDeltaK(nf,k,alphas,y,m_Q_ratio)
*
*     where:
*     - k   = 1: Hg, 2: Hq,
*     - as  = value of alphas / ( 4 * pi ),
*     - y   = Bjorken's variable.
*
************************************************************************
      function integrandsMatchingRes(y)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/grid.h"
      include "../commons/gridAlpha.h"
      include "../commons/wrapRes.h"
      include "../commons/ThresholdAlphaQCD.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      double precision y
**
*     Internal Variables
*
      double precision z,w_int,fR
      double precision xDeltaK,alphas,moQ
**
*     Output Variables
*
      double precision integrandsMatchingRes
*
*     Interpolant functions
*
      z = xg(igrid,wbeta) / y
*
      fR = w_int(inter_degree(igrid),walpha,z)
*
*     Contructing integrands
*
      alphas = 4d0 * pi * asthUp(wnf)
c      alphas = 4d0 * pi * asthDown(wnf)
      moQ = 1d0 / dsqrt(k2th(wnf))
      integrandsMatchingRes = xDeltaK(wnf,k,alphas,y,moQ) * fR
*
      return
      end

