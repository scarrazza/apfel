************************************************************************
*
*     integrandsDIS.f:
*
*     This functions return the integrands need to compute the DIS
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
*     Zero Mass coefficient functions
*
************************************************************************
      function integrandsDISzm(y)
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
      double precision C2G1A,C2NS1A,C2NS1B
      double precision CLG1A,CLNS1A
      double precision C3NS1A,C3NS1B
      double precision C2G2A,C2PS2A,C2NSP2A,C2NSM2A,C2NS2B
      double precision CLG2A,CLPS2A,CLNSP2A,CLNSM2A
      double precision C3NSP2A,C3NSM2A,C3NS2B
**
*     Output Variables
*
      double precision integrandsDISzm
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
               C2R(k,1) = C2G1A(y)
               C2S(k,1) = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               C2R(k,1) = 0d0
               C2S(k,1) = 0d0
*     Non-singlet-plus/minus
            elseif(k.eq.3.or.k.eq.4)then
               C2R(k,1) = C2NS1A(y)
               C2S(k,1) = C2NS1B(y)
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               CLR(k,1) = CLG1A(y)
               CLS(k,1) = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               CLR(k,1) = 0d0
               CLS(k,1) = 0d0
*     Non-singlet-plus/minus
            elseif(k.eq.3.or.k.eq.4)then
               CLR(k,1) = CLNS1A(y)
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
               C3R(k,1) = C3NS1A(y)
               C3S(k,1) = C3NS1B(y)
            endif
         endif
      endif
*
*     NNLO
*
      if(wipt.eq.2)then
*     C2
         if(sf.eq.1)then
*     Gluon
            if(k.eq.1)then
               C2R(k,2) = C2G2A(y,1)
               C2S(k,2) = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               C2R(k,2) = C2PS2A(y,1)
               C2S(k,2) = 0d0
*     Non-singlet-plus
            elseif(k.eq.3)then
               C2R(k,2) = C2NSP2A(y,wnf)
               C2S(k,2) = C2NS2B(y,wnf)
*     Non-singlet-minus
            elseif(k.eq.4)then
               C2R(k,2) = C2NSM2A(y,wnf)
               C2S(k,2) = C2NS2B(y,wnf)
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               CLR(k,2) = CLG2A(y,1)
               CLS(k,2) = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               CLR(k,2) = CLPS2A(y,1)
               CLS(k,2) = 0d0
*     Non-singlet-plus
            elseif(k.eq.3)then
               CLR(k,2) = CLNSP2A(y,wnf)
               CLS(k,2) = 0d0
*     Non-singlet-minus
            elseif(k.eq.4)then
               CLR(k,2) = CLNSM2A(y,wnf)
               CLS(k,2) = 0d0
            endif
*     C3
         elseif(sf.eq.3)then
*     Gluon
            if(k.eq.1)then
               C3R(k,2) = 0d0
               C3S(k,2) = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               C3R(k,2) = 0d0
               C3S(k,2) = 0d0
*     Non-singlet-plus
            elseif(k.eq.3)then
               C3R(k,2) = C3NSP2A(y,wnf)
               C3S(k,2) = C3NS2B(y,wnf)
*     Non-singlet-minus
            elseif(k.eq.4)then
               C3R(k,2) = C3NSM2A(y,wnf)
               C3S(k,2) = C3NS2B(y,wnf)
            endif
         endif
      endif
*
      if(sf.eq.1)then
         integrandsDISzm = C2R(k,wipt) * fR + C2S(k,wipt) * fS
c         integrandsDISzm = z * ( C2R(k,wipt) * fR + C2S(k,wipt) * fS ) / y
      elseif(sf.eq.2)then
         integrandsDISzm = CLR(k,wipt) * fR + CLS(k,wipt) * fS
c         integrandsDISzm = z * ( CLR(k,wipt) * fR + CLS(k,wipt) * fS ) / y
      elseif(sf.eq.3)then
         integrandsDISzm = C3R(k,wipt) * fR + C3S(k,wipt) * fS
c         integrandsDISzm = z * ( C3R(k,wipt) * fR + C3S(k,wipt) * fS ) / y
      endif
*
      return
      end
*
************************************************************************
*
*     Massive coefficient functions (Neutral Current)
*
************************************************************************
      function integrandsDISNCm(y)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/wrapDIS.h"
      include "../commons/coeffhqmellin.h"
**
*     Input Variables
*
      double precision y
**
*     Internal Variables
*
      double precision z,w_int,fR,fS,fL
      double precision C2R(3,2)
      double precision CLR(3,2)
      double precision C3R(3,2),C3S(3,2)
      double precision MassiveCF
      double precision C3NS1A,C3NS1B
      double precision C3NSP2A,C3NS2B
**
*     Output Variables
*
      double precision integrandsDISNCm
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
               C2R(k,1) = MassiveCF(1,wixi,y)
*     Pure-singlet
            elseif(k.eq.2)then
               C2R(k,1) = 0d0
*     Non-singlet
            elseif(k.eq.3)then
               C2R(k,1) = 0d0
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               CLR(k,1) = MassiveCF(2,wixi,y)
*     Pure-singlet
            elseif(k.eq.2)then
               CLR(k,1) = 0d0
*     Non-singlet
            elseif(k.eq.3)then
               CLR(k,1) = 0d0
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
*     Non-singlet
            elseif(k.eq.3)then
               C3R(k,1) = C3NS1A(y)
               C3S(k,1) = C3NS1B(y)
            endif
         endif
      endif
*
*     NNLO
*
      if(wipt.eq.2)then
*     C2
         if(sf.eq.1)then
*     Gluon
            if(k.eq.1)then
               C2R(k,2) = MassiveCF(3,wixi,y) 
     1                  + MassiveCF(9,wixi,y) * dlog(xigrid(wixi))
*     Pure-singlet
            elseif(k.eq.2)then
               C2R(k,2) = MassiveCF(4,wixi,y) 
     1                  + MassiveCF(10,wixi,y) * dlog(xigrid(wixi))
*     Non-singlet
            elseif(k.eq.3)then
               C2R(k,2) = MassiveCF(7,wixi,y) 
c     1                  + MassiveCF(13,wixi,y) * dlog(xigrid(wixi))
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               CLR(k,2) = MassiveCF(5,wixi,y) 
     1                  + MassiveCF(11,wixi,y) * dlog(xigrid(wixi))
*     Pure-singlet
            elseif(k.eq.2)then
               CLR(k,2) = MassiveCF(6,wixi,y) 
     1                  + MassiveCF(12,wixi,y) * dlog(xigrid(wixi))
*     Non-singlet
            elseif(k.eq.3)then
               CLR(k,2) = MassiveCF(8,wixi,y)
            endif
*     C3
         elseif(sf.eq.3)then
*     Gluon
            if(k.eq.1)then
               C3R(k,2) = 0d0
               C3S(k,2) = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               C3R(k,2) = 0d0
               C3S(k,2) = 0d0
*     Non-singlet
            elseif(k.eq.3)then
               C3R(k,2) = C3NSP2A(y,wnf)
               C3S(k,2) = C3NS2B(y,wnf)
            endif
         endif
      endif
*
      if(sf.eq.1)then
         integrandsDISNCm = C2R(k,wipt) * fR
c         integrandsDISNCm = z * C2R(k,wipt) * fR / y
      elseif(sf.eq.2)then
         integrandsDISNCm = CLR(k,wipt) * fR
c         integrandsDISNCm = z * CLR(k,wipt) / y
      elseif(sf.eq.3)then
         integrandsDISNCm = C3R(k,wipt) * fR + C3S(k,wipt) * fS
c         integrandsDISNCm = z * ( C3R(k,wipt) * fR + C3S(k,wipt) * fS ) / y
      endif
*
      return
      end
*
************************************************************************
*
*     Massive zero coefficient functions (Neutral Current)
*
*     wl =  1    2    3    4    5
*           A0   AQ   AQ2  AF   AQF
*
************************************************************************
      function integrandsDISNCm0(y)
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
      double precision C2R(3,2,5),C2S(3,2,5)
      double precision CLR(3,2,5)
      double precision C2G1AM0_A0,C2G1AM0_AQ,CLG1AM0_A0
      double precision C2G2AM0_A0,C2G2AM0_AQ,C2G2AM0_AQ2,C2G2AM0_AF
      double precision C2G2AM0_AQF
      double precision C2PS2AM0_A0,C2PS2AM0_AQ,C2PS2AM0_AQ2,C2PS2AM0_AF
      double precision C2PS2AM0_AQF
      double precision C2NS2AM0_A0,C2NS2AM0_AQ,C2NS2AM0_AQ2
      double precision C2NS2BM0_A0,C2NS2BM0_AQ,C2NS2BM0_AQ2
      double precision CLG2AM0_A0,CLG2AM0_AQ,CLG2AM0_AF
      double precision CLPS2AM0_A0,CLPS2AM0_AQ,CLPS2AM0_AF
      double precision CLNS2AM0_A0,CLNS2AM0_AQ
      double precision C3R(3,2),C3S(3,2)
      double precision C3NS1A,C3NS1B
      double precision C3NSP2A,C3NS2B
**
*     Output Variables
*
      double precision integrandsDISNCm0
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
               if(wl.eq.1)then
                  C2R(k,wipt,wl) = C2G1AM0_A0(y)
               elseif(wl.eq.2)then
                  C2R(k,wipt,wl) = C2G1AM0_AQ(y)
               else
                  C2R(k,wipt,wl) = 0d0
               endif
*     Pure-singlet
            elseif(k.eq.2)then
               C2R(k,wipt,wl) = 0d0
*     Non-singlet
            elseif(k.eq.3)then
               C2R(k,wipt,wl) = 0d0
            endif
            C2S(k,wipt,wl) = 0d0
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               if(wl.eq.1)then
                  CLR(k,wipt,wl) = CLG1AM0_A0(y)
               else
                  CLR(k,wipt,wl) = 0d0
               endif
*     Pure-singlet
            elseif(k.eq.2)then
               CLR(k,wipt,wl) = 0d0
*     Non-singlet
            elseif(k.eq.3)then
               CLR(k,wipt,wl) = 0d0
            endif
*     C3
         elseif(sf.eq.3)then
*     Gluon
            if(k.eq.1)then
               C3R(k,wipt) = 0d0
               C3S(k,wipt) = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               C3R(k,wipt) = 0d0
               C3S(k,wipt) = 0d0
*     Non-singlet
            elseif(k.eq.3)then
               C3R(k,wipt) = C3NS1A(y)
               C3S(k,wipt) = C3NS1B(y)
            endif
         endif
      endif
*
*     NNLO
*
      if(wipt.eq.2)then
*     C2
         if(sf.eq.1)then
*     Gluon
            if(k.eq.1)then
               if(wl.eq.1)then
                  C2R(k,wipt,wl) = C2G2AM0_A0(y)
               elseif(wl.eq.2)then
                  C2R(k,wipt,wl) = C2G2AM0_AQ(y)
               elseif(wl.eq.3)then
                  C2R(k,wipt,wl) = C2G2AM0_AQ2(y)
               elseif(wl.eq.4)then
                  C2R(k,wipt,wl) = C2G2AM0_AF(y)
               elseif(wl.eq.5)then
                  C2R(k,wipt,wl) = C2G2AM0_AQF(y)
               endif
               C2S(k,wipt,wl) = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               if(wl.eq.1)then
                  C2R(k,wipt,wl) = C2PS2AM0_A0(y)
               elseif(wl.eq.2)then
                  C2R(k,wipt,wl) = C2PS2AM0_AQ(y)
               elseif(wl.eq.3)then
                  C2R(k,wipt,wl) = C2PS2AM0_AQ2(y)
               elseif(wl.eq.4)then
                  C2R(k,wipt,wl) = C2PS2AM0_AF(y)
               elseif(wl.eq.5)then
                  C2R(k,wipt,wl) = C2PS2AM0_AQF(y)
               endif
               C2S(k,wipt,wl) = 0d0
*     Non-singlet
            elseif(k.eq.3)then
               if(wl.eq.1)then
                  C2R(k,wipt,wl) = C2NS2AM0_A0(y)
                  C2S(k,wipt,wl) = C2NS2BM0_A0(y)
               elseif(wl.eq.2)then
                  C2R(k,wipt,wl) = C2NS2AM0_AQ(y)
                  C2S(k,wipt,wl) = C2NS2BM0_AQ(y)
               elseif(wl.eq.3)then
                  C2R(k,wipt,wl) = C2NS2AM0_AQ2(y)
                  C2S(k,wipt,wl) = C2NS2BM0_AQ2(y)
               else
                  C2R(k,wipt,wl) = 0d0
                  C2S(k,wipt,wl) = 0d0
               endif
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               if(wl.eq.1)then
                  CLR(k,wipt,wl) = CLG2AM0_A0(y)
               elseif(wl.eq.2)then
                  CLR(k,wipt,wl) = CLG2AM0_AQ(y)
               elseif(wl.eq.4)then
                  CLR(k,wipt,wl) = CLG2AM0_AF(y)
               else
                  CLR(k,wipt,wl) = 0d0
               endif
*     Pure-singlet
            elseif(k.eq.2)then
               if(wl.eq.1)then
                  CLR(k,wipt,wl) = CLPS2AM0_A0(y)
               elseif(wl.eq.2)then
                  CLR(k,wipt,wl) = CLPS2AM0_AQ(y)
               elseif(wl.eq.4)then
                  CLR(k,wipt,wl) = CLPS2AM0_AF(y)
               else
                  CLR(k,wipt,wl) = 0d0
               endif
*     Non-singlet
            elseif(k.eq.3)then
               if(wl.eq.1)then
                  CLR(k,wipt,wl) = CLNS2AM0_A0(y)
               elseif(wl.eq.2)then
                  CLR(k,wipt,wl) = CLNS2AM0_AQ(y)
               else
                  CLR(k,wipt,wl) = 0d0
               endif
            endif
*     C3
         elseif(sf.eq.3)then
*     Gluon
            if(k.eq.1)then
               C3R(k,wipt) = 0d0
               C3S(k,wipt) = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               C3R(k,wipt) = 0d0
               C3S(k,wipt) = 0d0
*     Non-singlet
            elseif(k.eq.3)then
               C3R(k,wipt) = C3NSP2A(y,wnf)
               C3S(k,wipt) = C3NS2B(y,wnf)
            endif
         endif
      endif
*
      if(sf.eq.1)then
         integrandsDISNCm0 = C2R(k,wipt,wl) * fR + C2S(k,wipt,wl) * fS
      elseif(sf.eq.2)then
         integrandsDISNCm0 = CLR(k,wipt,wl) * fR
      elseif(sf.eq.3)then
         integrandsDISNCm0 = C3R(k,wipt) * fR + C3S(k,wipt) * fS
      endif
*
      return
      end
*
************************************************************************
*
*     Massive coefficient functions (Charged Current)
*
************************************************************************
      function integrandsDISCCm(y)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/wrapDIS.h"
      include "../commons/coeffhqmellin.h"
**
*     Input Variables
*
      double precision y
**
*     Internal Variables
*
      double precision z,w_int,fR,fS,fL
      double precision xi
      double precision C2R(3,2),C2S(3,2)
      double precision CLR(3,2),CLS(3,2)
      double precision C3R(3,2),C3S(3,2)
      double precision c2ns1cca,c2ns1ccb,c2g1cca
      double precision clns1cca,clns1ccb,clg1cca
      double precision c3ns1cca,c3ns1ccb,c3g1cca
**
*     Output Variables
*
      double precision integrandsDISCCm
*
      xi = xigrid(wixi)
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
               C2R(k,wipt) = c2g1cca(xi,y)
               C2S(k,wipt) = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               C2R(k,wipt) = 0d0
               C2S(k,wipt) = 0d0
*     Non-singlet
            elseif(k.eq.3)then
               C2R(k,wipt) = c2ns1cca(xi,y)
               C2S(k,wipt) = c2ns1ccb(xi,y)
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               CLR(k,wipt) = clg1cca(xi,y)
               CLS(k,wipt) = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               CLR(k,wipt) = 0d0
               CLS(k,wipt) = 0d0
*     Non-singlet
            elseif(k.eq.3)then
               CLR(k,wipt) = clns1cca(xi,y)
               CLS(k,wipt) = clns1ccb(xi,y)
            endif
*     C3
         elseif(sf.eq.3)then
*     Gluon
            if(k.eq.1)then
               C3R(k,wipt) = c3g1cca(xi,y)
               C3S(k,wipt) = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               C3R(k,wipt) = 0d0
               C3S(k,wipt) = 0d0
*     Non-singlet
            elseif(k.eq.3)then
               C3R(k,wipt) = c3ns1cca(xi,y)
               C3S(k,wipt) = c3ns1ccb(xi,y)
            endif
         endif
      endif
*
*     NNLO
*
      if(wipt.eq.2)then
*     C2
         C2R(k,wipt) = 0d0
         C2S(k,wipt) = 0d0
*     CL
         CLR(k,wipt) = 0d0
         CLS(k,wipt) = 0d0
*     C3
         C3R(k,wipt) = 0d0
         C3S(k,wipt) = 0d0
      endif
*
      if(sf.eq.1)then
         integrandsDISCCm = C2R(k,wipt) * fR + C2S(k,wipt) * fS
      elseif(sf.eq.2)then
         integrandsDISCCm = CLR(k,wipt) * fR + CLS(k,wipt) * fS
      elseif(sf.eq.3)then
         integrandsDISCCm = C3R(k,wipt) * fR + C3S(k,wipt) * fS
      endif
*
      return
      end
*
************************************************************************
*
*     Massive zero coefficient functions (Charged Current)
*
************************************************************************
      function integrandsDISCCm0(y)
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
      double precision z,w_int,fR
      double precision CG1ACCM0_AL
**
*     Output Variables
*
      double precision integrandsDISCCm0
*
*     Interpolant functions
*
      z = xg(igrid,wbeta) / y
*
      fR = w_int(inter_degree(igrid),walpha,z)
*
*     Contructing integrand
*
      integrandsDISCCm0 = CG1ACCM0_AL(y) * fR
*
      return
      end
