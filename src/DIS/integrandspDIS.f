************************************************************************
*
*     integrandspDIS.f:
*
*     This functions return the integrands need to compute the polarised DIS
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
*        1          F2 (g4)
*        2          FL (gL)
*        3          F3 (g1)
*
*     that are contained in the common block wrapDIS.h
*
************************************************************************
*
*     Zero Mass coefficient functions for polarised DIS
*
************************************************************************
      function integrandspDISzm(y)
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
      double precision C3G1PA, C3NS1PA, C3NS1PB
      double precision CLNS1PA
      double precision C2NS1PA, C2NS1PB
**
*     Output Variables
*
      double precision integrandspDISzm
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
*     Non-singlet-plus/minus
            if(k.eq.3.or.k.eq.4)then
               CR = C2NS1PA(y)
               CS = C2NS1PB(y)
            endif
*     CL
         elseif(sf.eq.2)then
*     Non-singlet-plus/minus
            if(k.eq.3.or.k.eq.4)then
               CR = CLNS1PA(y)
            endif
*     C3
         elseif(sf.eq.3)then
*     Gluon
            if(k.eq.1)then
               CR = C3G1PA(y)
*     Non-singlet-plus/minus
            elseif(k.eq.3.or.k.eq.4)then
               CR = C3NS1PA(y)
               CS = C3NS1PB(y)
            endif
         endif
C $$$ *
C $$$ *     NNLO
C $$$ *
C $$$     elseif(wipt.eq.2)then
C $$$ *     C2
C $$$        if(sf.eq.1)then
C $$$ *     Gluon
C $$$           if(k.eq.1)then
C $$$              CR = C2G2TA(y,1)
C $$$ *     Pure-singlet
C $$$           elseif(k.eq.2)then
C $$$              CR = C2PS2TA(y,1)
C $$$ *     Non-singlet-plus
C $$$           elseif(k.eq.3)then
C $$$              CR = C2NSP2TA(y,wnf)
C $$$              CS = C2NS2TB(y,wnf)
C $$$ *     Non-singlet-minus
C $$$           elseif(k.eq.4)then
C $$$              CR = C2NSP2TA(y,wnf)
C $$$              CS = C2NS2TB(y,wnf)
C $$$           endif
C $$$ *     CL
C $$$        elseif(sf.eq.2)then
C $$$ *     Gluon
C $$$           if(k.eq.1)then
C $$$              CR = CLG2TA(y,1)
C $$$ *     Pure-singlet
C $$$           elseif(k.eq.2)then
C $$$              CR = CLPS2TA(y,1)
C $$$ *     Non-singlet-plus
C $$$           elseif(k.eq.3)then
C $$$              CR = CLNSP2TA(y,wnf)
C $$$ *     Non-singlet-minus
C $$$           elseif(k.eq.4)then
C $$$              CR = CLNSP2TA(y,wnf)
C $$$           endif
C $$$ *     C3
C $$$        elseif(sf.eq.3)then
C $$$ *     Non-singlet-plus
C $$$           if(k.eq.3)then
C $$$              CR = C3NSP2TA(y,wnf)
C $$$              CS = C3NS2TB(y,wnf)
C $$$ *     Non-singlet-minus
C $$$           elseif(k.eq.4)then
C $$$              CR = C3NSP2TA(y,wnf)
C $$$              CS = C3NS2TB(y,wnf)
C $$$           endif
C $$$        endif
      endif
*
      integrandspDISzm = CR * fR + CS * fS
*
      return
      end
