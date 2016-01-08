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
      double precision CR,CS
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
      CR = 0d0
      CS = 0d0
*
*     NLO
*
      if(wipt.eq.1)then
*     C2
         if(sf.eq.1)then
*     Gluon
            if(k.eq.1)then
               CR = C2G1TA(y)
*     Non-singlet-plus/minus
            elseif(k.eq.3.or.k.eq.4)then
               CR = C2NS1TA(y)
               CS = C2NS1TB(y)
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               CR = CLG1TA(y)
*     Non-singlet-plus/minus
            elseif(k.eq.3.or.k.eq.4)then
               CR = CLNS1TA(y)
            endif
*     C3
         elseif(sf.eq.3)then
*     Non-singlet-plus/minus
            if(k.eq.3.or.k.eq.4)then
               CR = C3NS1TA(y)
               CS = C3NS1TB(y)
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
               CR = C2G2TA(y,1)
*     Pure-singlet
            elseif(k.eq.2)then
               CR = C2PS2TA(y,1)
*     Non-singlet-plus
            elseif(k.eq.3)then
               CR = C2NSP2TA(y,wnf)
               CS = C2NS2TB(y,wnf)
*     Non-singlet-minus
            elseif(k.eq.4)then
               CR = C2NSP2TA(y,wnf)
               CS = C2NS2TB(y,wnf)
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               CR = CLG2TA(y,1)
*     Pure-singlet
            elseif(k.eq.2)then
               CR = CLPS2TA(y,1)
*     Non-singlet-plus
            elseif(k.eq.3)then
               CR = CLNSP2TA(y,wnf)
*     Non-singlet-minus
            elseif(k.eq.4)then
               CR = CLNSP2TA(y,wnf)
            endif
*     C3
         elseif(sf.eq.3)then
*     Non-singlet-plus
            if(k.eq.3)then
               CR = C3NSP2TA(y,wnf)
               CS = C3NS2TB(y,wnf)
*     Non-singlet-minus
            elseif(k.eq.4)then
               CR = C3NSP2TA(y,wnf)
               CS = C3NS2TB(y,wnf)
            endif
         endif
      endif
*
      integrandsSIAzm = CR * fR + CS * fS
*
      return
      end
