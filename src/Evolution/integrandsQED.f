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
*        k  =  1   2   3   4   5   6   7   8   9   10  11
*              nsp nsm gg  gq  gD  qg  qq  qD  Dg  Dq  DD
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
*     NS up, NS down, Quark-Quark, Quark-Delta, Delta-Quark, Delta-Delta
         if(k.eq.1.or.k.eq.2.or.k.eq.7.or.
     1      k.eq.8.or.k.eq.10.or.k.eq.11)then
            PR = X0NSA(y) / CF
            PS = X0NSB(y) / CF
*     Quark-Gamma, Delta-Gamma
         elseif(k.eq.6.or.k.eq.9)then
            PR = X0QGA(y,1) / TR
            PS = 0d0
*     Gluon-Quark
         elseif(k.eq.4.or.k.eq.5)then
            PR = X0GQA(y) / CF
            PS = 0d0
*     Gamma-Gamma
         elseif(k.eq.3)then
            PR = 0d0
            PS = 0d0
         endif
      endif
*
*     NLO (i.e. O(alpha_s alpha))
*
      if(wipt.eq.1)then

      endif
*
*     NNLO (i.e. O(alpha^2))
*
      if(wipt.eq.2)then

      endif
*
      integrandsQED = PR * fR + PS * fS
*
      return
      end
