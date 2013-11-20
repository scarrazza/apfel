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
*     4) the particular plitting function denoted by k such that:
*
*        k  =  1   2   3   4   5   6   7
*              +   -   V   qq  qg  gq  gg
*
*     that are contained in the common block wrap.h
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
      double precision PR(7,0:2),PS(7,0:2)
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
*     Plus, Minus, Valence, Quark-Quark
      if(k.eq.1.or.k.eq.2.or.k.eq.3.or.k.eq.4)then
         PR(k,0) = X0NSA(y)
         PS(k,0) = X0NSB(y)
*     Quark-Gluon
      elseif(k.eq.5)then
         PR(k,0) = X0QGA(y,wnf)
         PS(k,0) = 0d0
*     Gluon-Quark
      elseif(k.eq.6)then
         PR(k,0) = X0GQA(y)
         PS(k,0) = 0d0
*     Gluon-Gluon
      elseif(k.eq.7)then
         PR(k,0) = X0GGA(y)
         PS(k,0) = X0GGB(y)
      endif
*
*     NLO
*
*     Plus
      if(wipt.ge.1)then
         if(k.eq.1)then
            PR(k,1) = X1NSPA(y,wnf)
            PS(k,1) = X1NSB(y)
*     Minus
         elseif(k.eq.2)then
            PR(k,1) = X1NSMA(y,wnf)
            PS(k,1) = X1NSB(y)
*     Valence
         elseif(k.eq.3)then
            PR(k,1) = X1NSMA(y,wnf)
            PS(k,1) = X1NSB(y)
*     Quark-Quark
         elseif(k.eq.4)then
            PR(k,1) = X1NSPA(y,wnf) + X1PSA(y,wnf)
            PS(k,1) = X1NSB(y)
*     Quark-Gluon
         elseif(k.eq.5)then
            PR(k,1) = X1QGA(y,wnf)
            PS(k,1) = 0d0
*     Gluon-Quark
         elseif(k.eq.6)then
            PR(k,1) = X1GQA(y,wnf)
            PS(k,1) = 0d0
*     Gluon-Gluon
         elseif(k.eq.7)then
            PR(k,1) = X1GGA(y,wnf)
            PS(k,1) = X1GGB(y)
         endif
      endif
*
*     NNLO
*
*     Plus
      if(wipt.ge.2)then
         if(k.eq.1)then
            PR(k,2) = P2NSPA(y,wnf)
            PS(k,2) = P2NSB(y,wnf)
*     Minus
         elseif(k.eq.2)then
            PR(k,2) = P2NSMA(y,wnf)
            PS(k,2) = P2NSB(y,wnf)
*     Valence
         elseif(k.eq.3)then
            PR(k,2) = P2NSMA(y,wnf) + P2NSSA(y,wnf)
            PS(k,2) = P2NSB(y,wnf)
*     Quark-Quark
         elseif(k.eq.4)then
            PR(k,2) = P2NSPA(y,wnf) + P2PSA(y,wnf)
            PS(k,2) = P2NSB(y,wnf)
*     Quark-Gluon
         elseif(k.eq.5)then
            PR(k,2) = P2QGA(y,wnf)
            PS(k,2) = 0d0
*     Gluon-Quark
         elseif(k.eq.6)then
            PR(k,2) = P2GQA(y,wnf)
            PS(k,2) = 0d0
*     Gluon-Gluon
         elseif(k.eq.7)then
            PR(k,2) = P2GGA(y,wnf)
            PS(k,2) = P2GGB(y,wnf)
         endif
      endif
*
      integrandsQCD = PR(k,wipt) * fR + PS(k,wipt) * fS
c      integrandsQCD = z * ( PR(k,wipt) * fR + PS(k,wipt) * fS ) / y
*
      return
      end
      
