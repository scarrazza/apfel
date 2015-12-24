************************************************************************
*
*     integrandsIC.f:
*
*     This functions return the integrands need to compute the IC
*     contributions to the structure functions.
*     In order to wrap it in a form that can be fed to the integration 
*     function DGAUSS, it is an explicit function only of the integration 
*     variable y while it depends on the other variables:
*
*     1) the grid indices walpha and wbeta,
*     2) the structure function index:
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
*     IC massive coefficient functions
*
************************************************************************
      function integrandsICm(y)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/wrapIC.h"
      include "../commons/wrapDIS.h"
**
*     Input Variables
*
      double precision y
**
*     Internal Variables
*
      double precision z,w_int,fS,fL
      double precision CR,CL
      double precision c21ICR,cL1ICR,c31ICR
      double precision one
      parameter(one=0.99999999d0)
**
*     Output Variables
*
      double precision integrandsICm
*
*     Interpolant functions
*
      lambda = 1d0 / xigrid(wixi)
      eta = 2d0 / ( 1d0 + dsqrt( 1d0 + 4d0 * lambda ) )
*
      integrandsICm = 0d0
      if(y.ge.eta) return
*
*     Interpolant functions
*
      z = xg(igrid,wbeta) / y / eta
*
      fR = w_int(inter_degree(igrid),walpha,z)
      fL = w_int(inter_degree(igrid),walpha,xg(igrid,wbeta)/eta)
*
*     Contructing integrands
*
*     C2
      if(sf.eq.1)then
         CR = c21ICR(y)
         CL = c21ICR(one)
*     CL
      elseif(sf.eq.2)then
         CR = cL1ICR(y)
         CL = cL1ICR(one)
*     C3
      elseif(sf.eq.3)then
         CR = c31ICR(y)
         CL = c31ICR(one)
      endif
*
      integrandsICm = ( CR * fR + CL * fL ) / ( 1d0 - y )
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
      double precision C2R,C2S
      double precision CLR
      double precision C3R,C3S
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
*     NLO
*
      if(wipt.eq.1)then
*     C2
         if(sf.eq.1)then
*     Gluon
            if(k.eq.1)then
               if(wl.eq.1)then
                  C2R = C2G1AM0_A0(y)
               elseif(wl.eq.2)then
                  C2R = C2G1AM0_AQ(y)
               else
                  C2R = 0d0
               endif
*     Pure-singlet
            elseif(k.eq.2)then
               C2R = 0d0
*     Non-singlet
            elseif(k.eq.3)then
               C2R = 0d0
            endif
            C2S = 0d0
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               if(wl.eq.1)then
                  CLR = CLG1AM0_A0(y)
               else
                  CLR = 0d0
               endif
*     Pure-singlet
            elseif(k.eq.2)then
               CLR = 0d0
*     Non-singlet
            elseif(k.eq.3)then
               CLR = 0d0
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
*     Non-singlet
            elseif(k.eq.3)then
               C3R = C3NS1A(y)
               C3S = C3NS1B(y)
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
                  C2R = C2G2AM0_A0(y)
               elseif(wl.eq.2)then
                  C2R = C2G2AM0_AQ(y)
               elseif(wl.eq.3)then
                  C2R = C2G2AM0_AQ2(y)
               elseif(wl.eq.4)then
                  C2R = C2G2AM0_AF(y)
               elseif(wl.eq.5)then
                  C2R = C2G2AM0_AQF(y)
               endif
               C2S = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               if(wl.eq.1)then
                  C2R = C2PS2AM0_A0(y)
               elseif(wl.eq.2)then
                  C2R = C2PS2AM0_AQ(y)
               elseif(wl.eq.3)then
                  C2R = C2PS2AM0_AQ2(y)
               elseif(wl.eq.4)then
                  C2R = C2PS2AM0_AF(y)
               elseif(wl.eq.5)then
                  C2R = C2PS2AM0_AQF(y)
               endif
               C2S = 0d0
*     Non-singlet
            elseif(k.eq.3)then
               if(wl.eq.1)then
                  C2R = C2NS2AM0_A0(y)
                  C2S = C2NS2BM0_A0(y)
               elseif(wl.eq.2)then
                  C2R = C2NS2AM0_AQ(y)
                  C2S = C2NS2BM0_AQ(y)
               elseif(wl.eq.3)then
                  C2R = C2NS2AM0_AQ2(y)
                  C2S = C2NS2BM0_AQ2(y)
               else
                  C2R = 0d0
                  C2S = 0d0
               endif
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               if(wl.eq.1)then
                  CLR = CLG2AM0_A0(y)
               elseif(wl.eq.2)then
                  CLR = CLG2AM0_AQ(y)
               elseif(wl.eq.4)then
                  CLR = CLG2AM0_AF(y)
               else
                  CLR = 0d0
               endif
*     Pure-singlet
            elseif(k.eq.2)then
               if(wl.eq.1)then
                  CLR = CLPS2AM0_A0(y)
               elseif(wl.eq.2)then
                  CLR = CLPS2AM0_AQ(y)
               elseif(wl.eq.4)then
                  CLR = CLPS2AM0_AF(y)
               else
                  CLR = 0d0
               endif
*     Non-singlet
            elseif(k.eq.3)then
               if(wl.eq.1)then
                  CLR = CLNS2AM0_A0(y)
               elseif(wl.eq.2)then
                  CLR = CLNS2AM0_AQ(y)
               else
                  CLR = 0d0
               endif
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
*     Non-singlet
            elseif(k.eq.3)then
               C3R = C3NSP2A(y,wnf)
               C3S = C3NS2B(y,wnf)
            endif
         endif
      endif
*
      if(sf.eq.1)then
         integrandsDISNCm0 = C2R * fR + C2S * fS
      elseif(sf.eq.2)then
         integrandsDISNCm0 = CLR * fR
      elseif(sf.eq.3)then
         integrandsDISNCm0 = C3R * fR + C3S * fS
      endif
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
      double precision C2R,C2S
      double precision CLR,CLS
      double precision C3R,C3S
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
*     NLO
*
      if(wipt.eq.1)then
*     C2
         if(sf.eq.1)then
*     Gluon
            if(k.eq.1)then
               C2R = c2g1cca(xi,ym)
               C2S = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               C2R = 0d0
               C2S = 0d0
*     Non-singlet
            elseif(k.eq.3)then
               C2R = c2ns1cca(xi,ym)
               C2S = c2ns1ccb(xi,ym)
            endif
*     CL
         elseif(sf.eq.2)then
*     Gluon
            if(k.eq.1)then
               CLR = clg1cca(xi,ym)
               CLS = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               CLR = 0d0
               CLS = 0d0
*     Non-singlet
            elseif(k.eq.3)then
               CLR = clns1cca(xi,ym)
               CLS = clns1ccb(xi,ym)
            endif
*     C3
         elseif(sf.eq.3)then
*     Gluon
            if(k.eq.1)then
               C3R = c3g1cca(xi,ym)
               C3S = 0d0
*     Pure-singlet
            elseif(k.eq.2)then
               C3R = 0d0
               C3S = 0d0
*     Non-singlet
            elseif(k.eq.3)then
               C3R = c3ns1cca(xi,ym)
               C3S = c3ns1ccb(xi,ym)
            endif
         endif
*
*     NNLO
*
      elseif(wipt.eq.2)then
*     C2
         C2R = 0d0
         C2S = 0d0
*     CL
         CLR = 0d0
         CLS = 0d0
*     C3
         C3R = 0d0
         C3S = 0d0
      endif
*
      if(sf.eq.1)then
         integrandsDISCCm = C2R * fR + C2S * fS
      elseif(sf.eq.2)then
         integrandsDISCCm = CLR * fR + CLS * fS
      elseif(sf.eq.3)then
         integrandsDISCCm = C3R * fR + C3S * fS
      endif
*
      integrandsDISCCm = integrandsDISCCm / lambda
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
