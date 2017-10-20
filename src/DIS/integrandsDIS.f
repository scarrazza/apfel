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
*     Zero Mass coefficient functions for DIS
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
      double precision CR,CS
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
               CR = C2G1A(y)
*     Non-singlet-plus/minus
            elseif(k.eq.3.or.k.eq.4)then
               CR = C2NS1A(y)
               CS = C2NS1B(y)
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               CR = CLG1A(y)
*     Non-singlet-plus/minus
            elseif(k.eq.3.or.k.eq.4)then
               CR = CLNS1A(y)
            endif
*     C3
         elseif(sf.eq.3)then
*     Non-singlet-plus/minus
            if(k.eq.3.or.k.eq.4)then
               CR = C3NS1A(y)
               CS = C3NS1B(y)
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
               CR = C2G2A(y,1)
*     Pure-singlet
            elseif(k.eq.2)then
               CR = C2PS2A(y,1)
*     Non-singlet-plus
            elseif(k.eq.3)then
               CR = C2NSP2A(y,wnf)
               CS = C2NS2B(y,wnf)
*     Non-singlet-minus
            elseif(k.eq.4)then
               CR = C2NSM2A(y,wnf)
               CS = C2NS2B(y,wnf)
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               CR = CLG2A(y,1)
*     Pure-singlet
            elseif(k.eq.2)then
               CR = CLPS2A(y,1)
*     Non-singlet-plus
            elseif(k.eq.3)then
               CR = CLNSP2A(y,wnf)
*     Non-singlet-minus
            elseif(k.eq.4)then
               CR = CLNSM2A(y,wnf)
            endif
*     C3
         elseif(sf.eq.3)then
*     Non-singlet-plus
            if(k.eq.3)then
               CR = C3NSP2A(y,wnf)
               CS = C3NS2B(y,wnf)
*     Non-singlet-minus
            elseif(k.eq.4)then
               CR = C3NSM2A(y,wnf)
               CS = C3NS2B(y,wnf)
            endif
         endif
      endif
*
      integrandsDISzm = CR * fR + CS * fS
*
      return
      end
*
************************************************************************
*
*     Massive coefficient functions (Neutral Current) for DIS
*
************************************************************************
      function integrandsDISNCm(y)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/wrapDIS.h"
      include "../commons/coeffhqmellin.h"
      include "../commons/kfacQ.h"
**
*     Input Variables
*
      double precision y
**
*     Internal Variables
*
      double precision z,w_int,fR,fS,fL
      double precision CR,CS
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
               CR = MassiveCF(1,wixi,y)
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               CR = MassiveCF(2,wixi,y)
            endif
*     C3
         elseif(sf.eq.3)then
*     Non-singlet
            if(k.eq.3)then
               CR = C3NS1A(y)
               CS = C3NS1B(y)
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
               CR = MassiveCF(3,wixi,y) 
     1            + MassiveCF(9,wixi,y)
     2            * dlog( xigrid(wixi) * kfacQ )
*     Pure-singlet
            elseif(k.eq.2)then
               CR = MassiveCF(4,wixi,y) 
     1            + MassiveCF(10,wixi,y)
     2            * dlog( xigrid(wixi) * kfacQ )
*     Non-singlet
            elseif(k.eq.3)then
               CR = MassiveCF(7,wixi,y) 
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               CR = MassiveCF(5,wixi,y) 
     1            + MassiveCF(11,wixi,y)
     2            * dlog( xigrid(wixi) * kfacQ )
*     Pure-singlet
            elseif(k.eq.2)then
               CR = MassiveCF(6,wixi,y) 
     1            + MassiveCF(12,wixi,y)
     2            * dlog( xigrid(wixi) * kfacQ )
*     Non-singlet
            elseif(k.eq.3)then
               CR = MassiveCF(8,wixi,y)
            endif
*     C3
         elseif(sf.eq.3)then
*     Non-singlet
            if(k.eq.3)then
               CR = C3NSP2A(y,wnf)
               CS = C3NS2B(y,wnf)
            endif
         endif
      endif
*
      integrandsDISNCm = CR * fR
      if(sf.eq.3) integrandsDISNCm = integrandsDISNCm + CS * fS
*
      return
      end
*
************************************************************************
*
*     Massive zero coefficient functions (Neutral Current) for DIS
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
      double precision CR,CS
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
               if(wl.eq.1)then
                  CR = C2G1AM0_A0(y)
               elseif(wl.eq.2)then
                  CR = C2G1AM0_AQ(y)
               endif
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               if(wl.eq.1)then
                  CR = CLG1AM0_A0(y)
               endif
            endif
*     C3
         elseif(sf.eq.3)then
*     Non-singlet
            if(k.eq.3)then
               CR = C3NS1A(y)
               CS = C3NS1B(y)
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
               if(wl.eq.1)then
                  CR = C2G2AM0_A0(y)
               elseif(wl.eq.2)then
                  CR = C2G2AM0_AQ(y)
               elseif(wl.eq.3)then
                  CR = C2G2AM0_AQ2(y)
               elseif(wl.eq.4)then
                  CR = C2G2AM0_AF(y)
               elseif(wl.eq.5)then
                  CR = C2G2AM0_AQF(y)
               endif
*     Pure-singlet
            elseif(k.eq.2)then
               if(wl.eq.1)then
                  CR = C2PS2AM0_A0(y)
               elseif(wl.eq.2)then
                  CR = C2PS2AM0_AQ(y)
               elseif(wl.eq.3)then
                  CR = C2PS2AM0_AQ2(y)
               elseif(wl.eq.4)then
                  CR = C2PS2AM0_AF(y)
               elseif(wl.eq.5)then
                  CR = C2PS2AM0_AQF(y)
               endif
*     Non-singlet
            elseif(k.eq.3)then
               if(wl.eq.1)then
                  CR = C2NS2AM0_A0(y)
                  CS = C2NS2BM0_A0(y)
               elseif(wl.eq.2)then
                  CR = C2NS2AM0_AQ(y)
                  CS = C2NS2BM0_AQ(y)
               elseif(wl.eq.3)then
                  CR = C2NS2AM0_AQ2(y)
                  CS = C2NS2BM0_AQ2(y)
               endif
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               if(wl.eq.1)then
                  CR = CLG2AM0_A0(y)
               elseif(wl.eq.2)then
                  CR = CLG2AM0_AQ(y)
               elseif(wl.eq.4)then
                  CR = CLG2AM0_AF(y)
               endif
*     Pure-singlet
            elseif(k.eq.2)then
               if(wl.eq.1)then
                  CR = CLPS2AM0_A0(y)
               elseif(wl.eq.2)then
                  CR = CLPS2AM0_AQ(y)
               elseif(wl.eq.4)then
                  CR = CLPS2AM0_AF(y)
               endif
*     Non-singlet
            elseif(k.eq.3)then
               if(wl.eq.1)then
                  CR = CLNS2AM0_A0(y)
               elseif(wl.eq.2)then
                  CR = CLNS2AM0_AQ(y)
               endif
            endif
*     C3
         elseif(sf.eq.3)then
            if(k.eq.3)then
               CR = C3NSP2A(y,wnf)
               CS = C3NS2B(y,wnf)
            endif
         endif
      endif
*
      integrandsDISNCm0 = CR * fR + CS * fS
*
      return
      end
*
************************************************************************
*
*     Massive coefficient functions (Charged Current) for DIS
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
      double precision xi,lambda,ym
      double precision CR,CS
      double precision c2ns1cca,c2ns1ccb,c2g1cca
      double precision clns1cca,clns1ccb,clg1cca
      double precision c3ns1cca,c3ns1ccb,c3g1cca
**
*     Output Variables
*
      double precision integrandsDISCCm
*
      xi = xigrid(wixi)
      lambda = xi / ( 1d0 + xi )
*
      integrandsDISCCm = 0d0
      if(y.ge.lambda) return
*
*     Interpolant functions
*
      z = xg(igrid,wbeta) / y
*
      fL = w_int(inter_degree(igrid),walpha,xg(igrid,wbeta)/lambda)
      fR = w_int(inter_degree(igrid),walpha,z)
      fS = fR - fL
*
*     rescaled variable
*
      ym = y / lambda
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
               CR = c2g1cca(xi,ym)
*     Non-singlet
            elseif(k.eq.3)then
               CR = c2ns1cca(xi,ym)
               CS = c2ns1ccb(xi,ym)
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               CR = clg1cca(xi,ym)
*     Non-singlet
            elseif(k.eq.3)then
               CR = clns1cca(xi,ym)
               CS = clns1ccb(xi,ym)
            endif
*     C3
         elseif(sf.eq.3)then
*     Gluon
            if(k.eq.1)then
               CR = c3g1cca(xi,ym)
*     Non-singlet
            elseif(k.eq.3)then
               CR = c3ns1cca(xi,ym)
               CS = c3ns1ccb(xi,ym)
            endif
         endif
*
*     no NNLO yet
*
      endif
*
      integrandsDISCCm = ( CR * fR + CS * fS ) / lambda
*
      return
      end
*
************************************************************************
*
*     Massive zero coefficient functions (Charged Current) for DIS
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
*
************************************************************************
*
*     Integrands needed for the computation of the Target Mass
*     Corrections.
*
************************************************************************
      function integrandsDISTMC(y)
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
      double precision w_int
**
*     Output Variables
*
      double precision integrandsDISTMC
*
*     Contructing integrand
*
      integrandsDISTMC = w_int(1,walpha,y)
*
      return
      end
*
************************************************************************
*
*     Integrands for the small-x resummed coefficient functions.
*     The following routine calls the function "xDeltaC2" and "xDeltaCL"
*     that are the Fortran wrapper of the c++ function provided by the
*     Bonvini's code HELL.
*     The input variables of the function "xDeltaC2" and "xDeltaCL" are:
*
*        xDeltaC2(nf,k,alphas,y)
*        xDeltaCL(nf,k,alphas,y)
*
*     where:
*     - k   = 1: g, 2: pure singlet,
*     - as  = value of alphas / ( 4 * pi ),
*     - y   = Bjorken's variable.
*
************************************************************************
      function integrandsDISzmRes(y)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/grid.h"
      include "../commons/gridAlpha.h"
      include "../commons/wrapResDIS.h"
**
*     Input Variables
*
      double precision y
**
*     Internal Variables
*
      double precision z,w_int,fR
      double precision xDeltaC2,xDeltaCL,alphas
**
*     Output Variables
*
      double precision integrandsDISzmRes
*
*     Interpolant functions
*
      z = xg(igrid,wbeta) / y
*
      fR = w_int(inter_degree(igrid),walpha,z)
*
*     Contructing integrands
*
      alphas = 4d0 * pi * ag(wtau)
*     C2
      if(sf.eq.1)then
         integrandsDISzmRes = xDeltaC2(nfg(wtau),k,alphas,y) * fR
*     CL
      elseif(sf.eq.2)then
         integrandsDISzmRes = xDeltaCL(nfg(wtau),k,alphas,y) * fR
      else
         integrandsDISzmRes = 0d0
      endif
*
      return
      end
*
************************************************************************
      function integrandsDISNCcharmRes(y)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/grid.h"
      include "../commons/gridAlpha.h"
      include "../commons/wrapResDIS.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      double precision y
**
*     Internal Variables
*
      double precision z,w_int,fR
      double precision xDeltaMC2,xDeltaMCL,alphas
      double precision mh,m_Q_ratio,eta
**
*     Output Variables
*
      double precision integrandsDISNCcharmRes
*
*     Heavy quark mass
*
      mh = dsqrt(m2ph(4))
      m_Q_ratio = mh / qag(wtau)
*
*     Return zero in the prohibited kinematic region
*
      integrandsDISNCcharmRes = 0d0
      eta = 1d0 / ( 1d0 + 4d0 * m_Q_ratio**2 )
      if(y.ge.eta) return
*
*     Interpolant functions
*
      z = xg(igrid,wbeta) / y
*
      fR = w_int(inter_degree(igrid),walpha,z)
*
*     Contructing integrands
*
      alphas = 4d0 * pi * ag(wtau)
*     C2
      if(sf.eq.1)then
         integrandsDISNCcharmRes =
     1        xDeltaMC2(nfg(wtau),k,alphas,y,m_Q_ratio) * fR
*     CL
      elseif(sf.eq.2)then
         integrandsDISNCcharmRes =
     1        xDeltaMCL(nfg(wtau),k,alphas,y,m_Q_ratio) * fR
      else
         integrandsDISNCcharmRes = 0d0
      endif
*
      return
      end
*
************************************************************************
      function integrandsDISNCcharm0Res(y)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/grid.h"
      include "../commons/gridAlpha.h"
      include "../commons/wrapResDIS.h"
      include "../commons/m2th.h"
**
*     Input Variables
*
      double precision y
**
*     Internal Variables
*
      double precision z,w_int,fR
      double precision xDeltaK,alphas
      double precision mh,m_Q_ratio,eta
**
*     Output Variables
*
      double precision integrandsDISNCcharm0Res
*
*     Heavy quark mass
*
      mh = dsqrt(m2ph(4))
      m_Q_ratio = mh / qag(wtau)
*
*     Return zero in the prohibited kinematic region
*
      integrandsDISNCcharm0Res = 0d0
      eta = 1d0 / ( 1d0 + 4d0 * m_Q_ratio**2 )
      if(y.ge.eta) return
*
*     Interpolant functions
*
      z = xg(igrid,wbeta) / y
*
      fR = w_int(inter_degree(igrid),walpha,z)
*
*     Contructing integrands
*
      alphas = 4d0 * pi * ag(wtau)
*     C2
      if(sf.eq.1)then
         integrandsDISNCcharm0Res =
     1        xDeltaK(nfg(wtau),k,alphas,y,m_Q_ratio) * fR
*     CL
c      elseif(sf.eq.2)then
c         integrandsDISNCcharm0Res =
c     1        xDeltaK(nfg(wtau),k,alphas,y,m_Q_ratio) * fR
      else
         integrandsDISNCcharm0Res = 0d0
      endif
*
      return
      end

