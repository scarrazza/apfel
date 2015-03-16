************************************************************************
*
*     integrandsSIA.f:
*
*     This functions return the integrands need to compute the SIA
*     structure functions.
*     In order to wrap it in a form that can be fed to the integration 
*     function DGAUSS, it is an explicit function only of the integration 
*     variable y while it depends on the other variables:
*
*     1) the perturbatibe order pt, 
*     2) the number of active flavours wnf,
*     3) the grid indices walpha and wbeta,
*     4) the particular plitting function denoted by k such that:
*
*        k      combination
*     --------------------------
*        1         gluon  
*        2      pure-singlet
*        3    non-singlet-plus
*        4    non-singlet-minus
*
*     5) the structure function index:
*
*        sf  Structure Function
*     --------------------------
*        1          F2
*        2          FL
*        3          F3
*
*     that are contained in the common block wrapDIS.h
*
************************************************************************
*
*     Zero Mass coefficient functions for SIA
*
************************************************************************
      function integrandsSIAzm(y)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/wrapDIS.h"
**
*     Input Variables
*
      double precision y
**
*     Internal Variables
*
      double precision z,w_int,fR,fS,fL
      double precision C2R(4,2),C2S(4,2)
      double precision CLR(4,2),CLS(4,2)
      double precision C3R(4,2),C3S(4,2)
      double precision C2G1TA,C2NS1TA,C2NS1TB
      double precision CLG1TA,CLNS1TA
      double precision C3NS1TA,C3NS1TB
c      double precision C2G2TA,C2PS2TA,C2NSP2TA,C2NSM2TA,C2NS2TB
c      double precision CLG2TA,CLPS2TA,CLNSP2TA,CLNSM2TA
c      double precision C3NSP2TA,C3NSM2TA,C3NS2TB
**
*     Output Variables
*
      double precision integrandsSIAzm
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
*     C2
         if(sf.eq.1)then
*     Gluon
            if(k.eq.1)then
               C2R(k,1) = C2G1TA(y)
               C2S(k,1) = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               C2R(k,1) = 0d0
               C2S(k,1) = 0d0
*     Non-singlet-plus/minus
            elseif(k.eq.3.or.k.eq.4)then
               C2R(k,1) = C2NS1TA(y)
               C2S(k,1) = C2NS1TB(y)
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               CLR(k,1) = CLG1TA(y)
               CLS(k,1) = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               CLR(k,1) = 0d0
               CLS(k,1) = 0d0
*     Non-singlet-plus/minus
            elseif(k.eq.3.or.k.eq.4)then
               CLR(k,1) = CLNS1TA(y)
               CLS(k,1) = 0d0
            endif
*     C3
         elseif(sf.eq.3)then
*     Gluon
            if(k.eq.1)then
               C3R(k,1) = 0d0
               C3S(k,1) = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               C3R(k,1) = 0d0
               C3S(k,1) = 0d0
*     Non-singlet-plus/minus
            elseif(k.eq.3.or.k.eq.4)then
               C3R(k,1) = C3NS1TA(y)
               C3S(k,1) = C3NS1TB(y)
            endif
         endif
      endif
c$$$*
c$$$*     NNLO
c$$$*
c$$$      if(wipt.eq.2)then
c$$$*     C2
c$$$         if(sf.eq.1)then
c$$$*     Gluon
c$$$            if(k.eq.1)then
c$$$               C2R(k,2) = C2G2TA(y,1)
c$$$               C2S(k,2) = 0d0
c$$$*     Pure-singlet
c$$$            elseif(k.eq.2)then
c$$$               C2R(k,2) = C2PS2TA(y,1)
c$$$               C2S(k,2) = 0d0
c$$$*     Non-singlet-plus
c$$$            elseif(k.eq.3)then
c$$$               C2R(k,2) = C2NSP2TA(y,wnf)
c$$$               C2S(k,2) = C2NS2TB(y,wnf)
c$$$*     Non-singlet-minus
c$$$            elseif(k.eq.4)then
c$$$               C2R(k,2) = C2NSM2TA(y,wnf)
c$$$               C2S(k,2) = C2NS2TB(y,wnf)
c$$$            endif
c$$$*     CL
c$$$         elseif(sf.eq.2)then
c$$$*     Gluon
c$$$            if(k.eq.1)then
c$$$               CLR(k,2) = CLG2TA(y,1)
c$$$               CLS(k,2) = 0d0
c$$$*     Pure-singlet
c$$$            elseif(k.eq.2)then
c$$$               CLR(k,2) = CLPS2TA(y,1)
c$$$               CLS(k,2) = 0d0
c$$$*     Non-singlet-plus
c$$$            elseif(k.eq.3)then
c$$$               CLR(k,2) = CLNSP2TA(y,wnf)
c$$$               CLS(k,2) = 0d0
c$$$*     Non-singlet-minus
c$$$            elseif(k.eq.4)then
c$$$               CLR(k,2) = CLNSM2TA(y,wnf)
c$$$               CLS(k,2) = 0d0
c$$$            endif
c$$$*     C3
c$$$         elseif(sf.eq.3)then
c$$$*     Gluon
c$$$            if(k.eq.1)then
c$$$               C3R(k,2) = 0d0
c$$$               C3S(k,2) = 0d0
c$$$*     Pure-singlet
c$$$            elseif(k.eq.2)then
c$$$               C3R(k,2) = 0d0
c$$$               C3S(k,2) = 0d0
c$$$*     Non-singlet-plus
c$$$            elseif(k.eq.3)then
c$$$               C3R(k,2) = C3NSP2TA(y,wnf)
c$$$               C3S(k,2) = C3NS2TB(y,wnf)
c$$$*     Non-singlet-minus
c$$$            elseif(k.eq.4)then
c$$$               C3R(k,2) = C3NSM2TA(y,wnf)
c$$$               C3S(k,2) = C3NS2TB(y,wnf)
c$$$            endif
c$$$         endif
c$$$      endif
*
      if(sf.eq.1)then
         integrandsSIAzm = C2R(k,wipt) * fR + C2S(k,wipt) * fS
c         integrandsSIAzm = z * ( C2R(k,wipt) * fR + C2S(k,wipt) * fS ) / y
      elseif(sf.eq.2)then
         integrandsSIAzm = CLR(k,wipt) * fR + CLS(k,wipt) * fS
c         integrandsSIAzm = z * ( CLR(k,wipt) * fR + CLS(k,wipt) * fS ) / y
      elseif(sf.eq.3)then
         integrandsSIAzm = C3R(k,wipt) * fR + C3S(k,wipt) * fS
c         integrandsSIAzm = z * ( C3R(k,wipt) * fR + C3S(k,wipt) * fS ) / y
      endif
*
      return
      end
