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
      double precision C2R,C2S
      double precision CLR,CLS
      double precision C3R,C3S
      double precision C2G1TA,C2NS1TA,C2NS1TB
      double precision CLG1TA,CLNS1TA
      double precision C3NS1TA,C3NS1TB
      double precision C2G2TA,C2PS2TA,C2NSP2TA,C2NS2TB
      double precision CLG2TA,CLPS2TA,CLNSP2TA
      double precision C3NSP2TA,C3NS2TB
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
               C2R = C2G1TA(y)
               C2S = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               C2R = 0d0
               C2S = 0d0
*     Non-singlet-plus/minus
            elseif(k.eq.3.or.k.eq.4)then
               C2R = C2NS1TA(y)
               C2S = C2NS1TB(y)
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               CLR = CLG1TA(y)
               CLS = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               CLR = 0d0
               CLS = 0d0
*     Non-singlet-plus/minus
            elseif(k.eq.3.or.k.eq.4)then
               CLR = CLNS1TA(y)
               CLS = 0d0
            endif
*     C3
         elseif(sf.eq.3)then
*     Gluon
            if(k.eq.1)then
               C3R = 0d0
               C3S = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               C3R = 0d0
               C3S = 0d0
*     Non-singlet-plus/minus
            elseif(k.eq.3.or.k.eq.4)then
               C3R = C3NS1TA(y)
               C3S = C3NS1TB(y)
            endif
         endif
*
*     NNLO
*
      elseif(wipt.eq.2)then
*     C2
         if(sf.eq.1)then
*     Gluon
            if(k.eq.1)then
               C2R = C2G2TA(y,1)
               C2S = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               C2R = C2PS2TA(y,1)
               C2S = 0d0
*     Non-singlet-plus
            elseif(k.eq.3)then
               C2R = C2NSP2TA(y,wnf)
               C2S = C2NS2TB(y,wnf)
*     Non-singlet-minus
            elseif(k.eq.4)then
               C2R = C2NSP2TA(y,wnf)
               C2S = C2NS2TB(y,wnf)
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               CLR = CLG2TA(y,1)
               CLS = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               CLR = CLPS2TA(y,1)
               CLS = 0d0
*     Non-singlet-plus
            elseif(k.eq.3)then
               CLR = CLNSP2TA(y,wnf)
               CLS = 0d0
*     Non-singlet-minus
            elseif(k.eq.4)then
               CLR = CLNSP2TA(y,wnf)
               CLS = 0d0
            endif
*     C3
         elseif(sf.eq.3)then
*     Gluon
            if(k.eq.1)then
               C3R = 0d0
               C3S = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               C3R = 0d0
               C3S = 0d0
*     Non-singlet-plus
            elseif(k.eq.3)then
               C3R = C3NSP2TA(y,wnf)
               C3S = C3NS2TB(y,wnf)
*     Non-singlet-minus
            elseif(k.eq.4)then
               C3R = C3NSP2TA(y,wnf)
               C3S = C3NS2TB(y,wnf)
            endif
         endif
      endif
*
      if(sf.eq.1)then
         integrandsSIAzm = C2R * fR + C2S * fS
      elseif(sf.eq.2)then
         integrandsSIAzm = CLR * fR + CLS * fS
      elseif(sf.eq.3)then
         integrandsSIAzm = C3R * fR + C3S * fS
      endif
*
      return
      end
