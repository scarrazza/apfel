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
*        1          g1
*        2          gL
*        3          g4
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
      double precision G1G1PA, G1NS1PA, G1NS1PB
      double precision GLNS1PA
      double precision G4NS1PA, G4NS1PB
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
*     g4
         if(sf.eq.3)then
*     Non-singlet-plus/minus
            if(k.eq.3.or.k.eq.4)then
               CR = G4NS1PA(y)
               CS = G4NS1PB(y)
            endif
*     gL
         elseif(sf.eq.2)then
*     Non-singlet-plus/minus
            if(k.eq.3.or.k.eq.4)then
               CR = GLNS1PA(y)
            endif
*     g1
         elseif(sf.eq.1)then
*     Gluon
            if(k.eq.1)then
               CR = G1G1PA(y)
*     Non-singlet-plus/minus
            elseif(k.eq.3.or.k.eq.4)then
               CR = G1NS1PA(y)
               CS = G1NS1PB(y)
            endif
         endif
      endif
*
      integrandspDISzm = CR * fR + CS * fS
*
      return
      end
