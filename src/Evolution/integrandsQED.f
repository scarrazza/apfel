************************************************************************
*
*     integrandsQED.f:
*
*     This function returns the integrands need to compute the evolution 
*     operators in QED.
*     In order to wrap it in a form that can be fed to the integration 
*     function DGAUSS, it is an explicit function only of the integration 
*     variable y while it depends on the other variables:
*
*     1) the perturbatibe order pt, 
*     2) the grid indices walpha and wbeta,
*     3) the particular splitting function denoted by k such that:
*
*     k  =  1    2    3    4
*           nspu nsmu nspd nsmd
*
*           5    6    7    8
*           gg   ggm  gS   gD
*
*           9    10   11   12
*           gmg  gmgm gmS  gmD
*
*           13   14   15   16
*           Sg   Sgm  SS   SD
*
*           17   18   19   20
*           Dg   Dgm  DS   DD
*
*           21   22
*           VV   VDV (DVDV = VV, DVV = VDV)
*
*           23   24   25
*           LL   gmL  Lgm
*
*     that are contained in the common block wrap.h.
* 
************************************************************************
      function integrandsQED(y)
*
      implicit none
*
      include "../commons/ColorFactors.h"
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
      double precision X0NSA,X0NSB,X0QGA,X0GQA
      double precision X1NSPA_ASA,X1NSMA_ASA,X1QGAMA_ASA,X1GAMQA_ASA
      double precision X1GGAMA_ASA
      double precision REM_X1NSA_AA,REM_X1NSB_AA,REM_X1GAMQA_AA,X1PSA_AA
**
*     Output Variables
*
      double precision integrandsQED
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
*     Contructing integrands
*
*     LO
*
      if(wipt.eq.0)then
*     NS+ up, NS+ down, NS- up, NS- down, Quark-Quark, Quark-Delta, Delta-Quark, Delta-Delta, Valence-Valence, Valence-D_V
         if(k.eq.1.or.k.eq.2.or.k.eq.3.or.k.eq.4.or.
     1      k.eq.15.or.k.eq.16.or.k.eq.19.or.k.eq.20.or.
     2      k.eq.21.or.k.eq.22)then
            PR = X0NSA(y) / CF
            PS = X0NSB(y) / CF
*     Quark-Gamma, Delta-Gamma
         elseif(k.eq.14.or.k.eq.18)then
            PR = X0QGA(y,1) / TR
            PS = 0d0
*     Gamma-Quark, Gamma-Delta
         elseif(k.eq.11.or.k.eq.12)then
            PR = X0GQA(y) / CF
            PS = 0d0
*     all others are null
         else
            PR = 0d0
            PS = 0d0
         endif
      endif
*
*     NLO (i.e. O(alpha_s alpha))
*
      if(wipt.eq.1)then
*     NS+ up, NS+ down, Quark-Quark, Quark-Delta, Delta-Quark, Delta-Delta
         if(k.eq.1.or.k.eq.2.or.
     1      k.eq.15.or.k.eq.16.or.k.eq.19.or.k.eq.20)then
            PR = X1NSPA_ASA(y)
            PS = 0d0
*     NS- up, NS- down, Valence-Valence, Valence-D_V
         elseif(k.eq.3.or.k.eq.4.or.
     1          k.eq.21.or.k.eq.22)then
            PR = X1NSMA_ASA(y)
            PS = 0d0
*     Quark-Gluon, Delta-Gluon, Quark-Gamma, Delta-Gamma
         elseif(k.eq.13.or.k.eq.14.or.k.eq.17.or.k.eq.18)then
            PR = X1QGAMA_ASA(y)
            PS = 0d0
*     Gluon-Quark, Gluon-Delta, Gamma-Quark, Gamma-Delta
         elseif(k.eq.7.or.k.eq.8.or.k.eq.11.or.k.eq.12)then
            PR = X1GAMQA_ASA(y)
            PS = 0d0
*     Gluon-Gamma, Gamma-Gluon
         elseif(k.eq.6.or.k.eq.9)then
            PR = X1GGAMA_ASA(y)
            PS = 0d0
*     all others are null
         else
            PR = 0d0
            PS = 0d0
         endif
      endif
*
*     NNLO (i.e. O(alpha^2))
*
      if(wipt.eq.2)then
*     NS+ up, NS+ down, NS- up, NS- down
         if(k.eq.1.or.k.eq.2.or.k.eq.3.or.k.eq.4)then
            PR = REM_X1NSA_AA(y)
            PS = REM_X1NSB_AA(y)
*     Gamma-Quark, Gamma-Delta
         elseif(k.eq.11.or.k.eq.12)then
            PR = REM_X1GAMQA_AA(y)
            PS = 0d0
         elseif(k.eq.15.or.k.eq.16.or.k.eq.19.or.k.eq.20)then
            PR = X1PSA_AA(y)
            PS = 0d0
*     all others are null
         else
            PR = 0d0
            PS = 0d0
         endif
      endif
*
      integrandsQED = PR * fR + PS * fS
*
      return
      end
