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
*        k =    1            2            3
*             gluon     non-singlet  pure-singlet
*
*     5) the structure function index:
*
*        sf = 1   2   3
*             F2  FL  F3
*
*     that are contained in the common block wrapDIS.h
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
      double precision C2R(3,2),C2S(3,2)
      double precision CLR(3,2),CLS(3,2)
      double precision C3R(3,2),C3S(3,2)
      double precision C2G1A,C2NS1A,C2NS1B
      double precision CLG1A,CLNS1A
      double precision C3NS1A,C3NS1B
      double precision C2G2A,C2PS2A,C2NSP2A,C2NS2B
      double precision CLG2A,CLPS2A,CLNSP2A
      double precision C3NSP2A,C3NS2B
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
      if(wipt.ge.1)then
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
*     Non-singlet
            elseif(k.eq.3)then
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
*     Non-singlet
            elseif(k.eq.3)then
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
      if(wipt.ge.2)then
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
*     Non-singlet
            elseif(k.eq.3)then
               C2R(k,2) = C2NSP2A(y,wnf)
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
*     Non-singlet
            elseif(k.eq.3)then
               CLR(k,2) = CLNSP2A(y,wnf)
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
*     Non-singlet
            elseif(k.eq.3)then
               C3R(k,2) = C3NSP2A(y,wnf)
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
      function integrandsDISm(y)
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
      double precision integrandsDISm
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
      if(wipt.ge.1)then
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
      if(wipt.ge.2)then
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
               CLR(k,2) = MassiveCF(7,wixi,y)
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
         integrandsDISm = C2R(k,wipt) * fR
c         integrandsDISm = z * C2R(k,wipt) * fR / y
      elseif(sf.eq.2)then
         integrandsDISm = CLR(k,wipt) * fR
c         integrandsDISm = z * CLR(k,wipt) / y
      elseif(sf.eq.3)then
         integrandsDISm = C3R(k,wipt) * fR + C3S(k,wipt) * fS
c         integrandsDISm = z * ( C3R(k,wipt) * fR + C3S(k,wipt) * fS ) / y
      endif
*
      return
      end


