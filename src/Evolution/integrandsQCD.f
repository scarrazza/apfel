************************************************************************
*
*     integrandsQCD.f:
*
*     This function returns the integrands need to compute the evolution 
*     operators in QCD.
*     In order to wrap it in a form that can be fed to the integration 
*     function DGAUSS, it is an explicit function only of the integration 
*     variable y while it depends on the other variables:
*
*     1) the perturbatibe order pt, 
*     2) the number of active flavours wnf,
*     3) the grid indices walpha and wbeta,
*     4) the particular splitting function denoted by k such that:
*
*        k  =  1   2   3   4   5   6   7
*              +   -   V   qq  qg  gq  gg
*
*     that are contained in the common block wrap.h.
*
************************************************************************
      function integrandsQCD(y)
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
      double precision PR,PS
      double precision X0NSA,X0NSB,X0QGA,X0GQA,X0GGA,X0GGB
      double precision X1NSPA,X1NSB,X1NSMA,X1PSA,X1QGA,X1GQA,X1GGA,X1GGB
      double precision P2NSPA,P2NSB,P2NSMA,P2NSSA,P2PSA,P2QGA,P2GQA
      double precision P2GGA,P2GGB
**
*     Output Variables
*
      double precision integrandsQCD
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
*     LO
*
      if(wipt.eq.0)then
*     Plus, Minus, Valence, Quark-Quark
         if(k.eq.1.or.k.eq.2.or.k.eq.3.or.k.eq.4)then
            PR = X0NSA(y)
            PS = X0NSB(y)
*     Quark-Gluon
         elseif(k.eq.5)then
            PR = X0QGA(y,wnf)
            PS = 0d0
*     Gluon-Quark
         elseif(k.eq.6)then
            PR = X0GQA(y)
            PS = 0d0
*     Gluon-Gluon
         elseif(k.eq.7)then
            PR = X0GGA(y)
            PS = X0GGB(y)
         endif
*
*     NLO
*
      elseif(wipt.eq.1)then
*     Plus
         if(k.eq.1)then
            PR = X1NSPA(y,wnf)
            PS = X1NSB(y)
*     Minus
         elseif(k.eq.2)then
            PR = X1NSMA(y,wnf)
            PS = X1NSB(y)
*     Valence
         elseif(k.eq.3)then
            PR = X1NSMA(y,wnf)
            PS = X1NSB(y)
*     Quark-Quark
         elseif(k.eq.4)then
            PR = X1NSPA(y,wnf) + X1PSA(y,wnf)
            PS = X1NSB(y)
*     Quark-Gluon
         elseif(k.eq.5)then
            PR = X1QGA(y,wnf)
            PS = 0d0
*     Gluon-Quark
         elseif(k.eq.6)then
            PR = X1GQA(y,wnf)
            PS = 0d0
*     Gluon-Gluon
         elseif(k.eq.7)then
            PR = X1GGA(y,wnf)
            PS = X1GGB(y)
         endif
*
*     NNLO
*
      elseif(wipt.eq.2)then
*     Plus
         if(k.eq.1)then
            PR = P2NSPA(y,wnf)
            PS = P2NSB(y,wnf)
*     Minus
         elseif(k.eq.2)then
            PR = P2NSMA(y,wnf)
            PS = P2NSB(y,wnf)
*     Valence
         elseif(k.eq.3)then
            PR = P2NSMA(y,wnf) + P2NSSA(y,wnf)
            PS = P2NSB(y,wnf)
*     Quark-Quark
         elseif(k.eq.4)then
            PR = P2NSPA(y,wnf) + P2PSA(y,wnf)
            PS = P2NSB(y,wnf)
*     Quark-Gluon
         elseif(k.eq.5)then
            PR = P2QGA(y,wnf)
            PS = 0d0
*     Gluon-Quark
         elseif(k.eq.6)then
            PR = P2GQA(y,wnf)
            PS = 0d0
*     Gluon-Gluon
         elseif(k.eq.7)then
            PR = P2GGA(y,wnf)
            PS = P2GGB(y,wnf)
         endif
      endif
*
      integrandsQCD = PR * fR + PS * fS
*
      return
      end
*
************************************************************************
*
*     Integrands for the time-like evolution.
*
************************************************************************
      function integrandsQCDT(y)
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
      double precision PR,PS
      double precision X0NSA,X0NSB,X0QGA,X0GQA,X0GGA,X0GGB
      double precision X1NSPTA,X1NSTB,X1NSMTA,X1PSTA,X1QGTA,X1GQTA
      double precision X1GGTA,X1GGTB
      double precision P2NSPTA,P2NSTB,P2NSMTA,P2NSSTA
      double precision P2GGTA,P2GGTB,P2PSTA,P2QGTA,P2GQTA
**
*     Output Variables
*
      double precision integrandsQCDT
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
*     LO
*
      if(wipt.eq.0)then
*     Plus, Minus, Valence, Quark-Quark
         if(k.eq.1.or.k.eq.2.or.k.eq.3.or.k.eq.4)then
            PR = X0NSA(y)
            PS = X0NSB(y)
*     Quark-Gluon
         elseif(k.eq.5)then
            PR = X0QGA(y,wnf)
            PS = 0d0
*     Gluon-Quark
         elseif(k.eq.6)then
            PR = X0GQA(y)
            PS = 0d0
*     Gluon-Gluon
         elseif(k.eq.7)then
            PR = X0GGA(y)
            PS = X0GGB(y)
         endif
*
*     NLO
*
      elseif(wipt.eq.1)then
*     Plus
         if(k.eq.1)then
            PR = X1NSPTA(y,wnf)
            PS = X1NSTB(y)
*     Minus
         elseif(k.eq.2)then
            PR = X1NSMTA(y,wnf)
            PS = X1NSTB(y)
*     Valence
         elseif(k.eq.3)then
            PR = X1NSMTA(y,wnf)
            PS = X1NSTB(y)
*     Quark-Quark
         elseif(k.eq.4)then
            PR = X1NSPTA(y,wnf) + X1PSTA(y,wnf)
            PS = X1NSTB(y)
*     Quark-Gluon
         elseif(k.eq.5)then
            PR = X1QGTA(y,wnf)
            PS = 0d0
*     Gluon-Quark
         elseif(k.eq.6)then
            PR = X1GQTA(y,wnf)
            PS = 0d0
*     Gluon-Gluon
         elseif(k.eq.7)then
            PR = X1GGTA(y,wnf)
            PS = X1GGTB(y)
         endif
*
*     NNLO
*
      elseif(wipt.eq.2)then
*     Plus
         if(k.eq.1)then
            PR = P2NSPTA(y,wnf)
            PS = P2NSTB(y,wnf)
*     Minus
         elseif(k.eq.2)then
            PR = P2NSMTA(y,wnf)
            PS = P2NSTB(y,wnf)
*     Valence
         elseif(k.eq.3)then
            PR = P2NSMTA(y,wnf) + P2NSSTA(y,wnf)
            PS = P2NSTB(y,wnf)
*     Quark-Quark
         elseif(k.eq.4)then
            PR = P2NSPTA(y,wnf) + P2PSTA(y,wnf)
            PS = P2NSTB(y,wnf)
*     Quark-Gluon
         elseif(k.eq.5)then
            PR = P2QGTA(y,wnf)
            PS = 0d0
*     Gluon-Quark
         elseif(k.eq.6)then
            PR = P2GQTA(y,wnf)
            PS = 0d0
*     Gluon-Gluon
         elseif(k.eq.7)then
            PR = P2GGTA(y,wnf)
            PS = P2GGTB(y,wnf)
         endif
      endif
*
      integrandsQCDT = PR * fR + PS * fS
*
      return
      end
*
************************************************************************
*
*     Integrands for the small-x resummed evolution.
*     The following routine calls the function "xDeltaP" that is the
*     Fortran wrapper of the c++ function provided by the Bonvini's
*     code HELL.
*     The input variables of the function "xDeltaP" are:
*
*        xDeltaP(nf,k,alphas,y)
*
*     where:
*     - k   = 4: qq, 5: qg, 6: gq, 7: gg (splitting function),
*     - as  = value of alphas / ( 4 * pi ),
*     - y   = Bjorken's variable.
*
************************************************************************
      function integrandsQCDRes(y)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/grid.h"
      include "../commons/gridAlpha.h"
      include "../commons/wrapRes.h"
**
*     Input Variables
*
      double precision y
**
*     Internal Variables
*
      double precision z,w_int,fR
      double precision xDeltaP,alphas
**
*     Output Variables
*
      double precision integrandsQCDRes
*
*     Interpolant functions
*
      z = xg(igrid,wbeta) / y
*
      fR = w_int(inter_degree(igrid),walpha,z)
*
*     Contructing integrands
*
      alphas = 4d0 * pi * ag(wtau)
      integrandsQCDRes = xDeltaP(nfg(wtau),k,alphas,y) * fR
*
      return
      end
*
************************************************************************
*
*     Integrands for the polarized evolution.
*
************************************************************************
      function integrandsQCDPol(y)
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
      double precision PR,PS
      double precision X0NSPA,X0NSPB,X0QGPA,X0GQPA,X0GGPA,X0GGPB
      double precision X1NSPPA,X1NSPB,X1NSMPA,X1PSPA,X1QGPA,X1GQPA
      double precision X1GGPA,X1GGPB
      double precision P2NSPPA,P2NSPB,P2NSMPA,P2NSSPA,P2PSPA,P2QGPA
      double precision P2GQPA,P2GGPA,P2GGPB
**
*     Output Variables
*
      double precision integrandsQCDPol
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
*     LO
*
      if(wipt.eq.0)then
*     Plus, Minus, Valence, Quark-Quark
         if(k.eq.1.or.k.eq.2.or.k.eq.3.or.k.eq.4)then
            PR = X0NSPA(y)
            PS = X0NSPB(y)
*     Quark-Gluon
         elseif(k.eq.5)then
            PR = X0QGPA(y,wnf)
            PS = 0d0
*     Gluon-Quark
         elseif(k.eq.6)then
            PR = X0GQPA(y)
            PS = 0d0
*     Gluon-Gluon
         elseif(k.eq.7)then
            PR = X0GGPA(y)
            PS = X0GGPB(y)
         endif
*
*     NLO
*
      elseif(wipt.eq.1)then
*     Plus
         if(k.eq.1)then
            PR = X1NSPPA(y,wnf)
            PS = X1NSPB(y)
*     Minus
         elseif(k.eq.2)then
            PR = X1NSMPA(y,wnf)
            PS = X1NSPB(y)
*     Valence
         elseif(k.eq.3)then
            PR = X1NSMPA(y,wnf)
            PS = X1NSPB(y)
*     Quark-Quark
         elseif(k.eq.4)then
            PR = X1NSPPA(y,wnf) + X1PSPA(y,wnf)
            PS = X1NSPB(y)
*     Quark-Gluon
         elseif(k.eq.5)then
            PR = X1QGPA(y,wnf)
            PS = 0d0
*     Gluon-Quark
         elseif(k.eq.6)then
            PR = X1GQPA(y,wnf)
            PS = 0d0
*     Gluon-Gluon
         elseif(k.eq.7)then
            PR = X1GGPA(y,wnf)
            PS = X1GGPB(y)
         endif
*
*     NNLO
*
      elseif(wipt.eq.2)then
*     Plus
         if(k.eq.1)then
            PR = P2NSPPA(y,wnf)
            PS = P2NSPB(y,wnf)
*     Minus
         elseif(k.eq.2)then
            PR = P2NSMPA(y,wnf)
            PS = P2NSPB(y,wnf)
*     Valence
         elseif(k.eq.3)then
            PR = P2NSMPA(y,wnf) + P2NSSPA(y,wnf)
            PS = P2NSPB(y,wnf)
*     Quark-Quark
         elseif(k.eq.4)then
            PR = P2NSPPA(y,wnf) + P2PSPA(y,wnf)
            PS = P2NSPB(y,wnf)
*     Quark-Gluon
         elseif(k.eq.5)then
            PR = P2QGPA(y,wnf)
            PS = 0d0
*     Gluon-Quark
         elseif(k.eq.6)then
            PR = P2GQPA(y,wnf)
            PS = 0d0
*     Gluon-Gluon
         elseif(k.eq.7)then
            PR = P2GGPA(y,wnf)
            PS = P2GGPB(y,wnf)
         endif
      endif
*
      integrandsQCDPol = PR * fR + PS * fS
*
      return
      end
