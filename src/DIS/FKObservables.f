************************************************************************
*
*     FKObservables.f:
*
*     This function returns simulates the FKgenerator of the NNPDF code.
*
************************************************************************
      function FKObservables(x,Q,y)
*
      implicit none
*
      include "../commons/consts.h"
      include "../commons/FKObservable.h"
      include "../commons/ipt.h"
      include "../commons/eftcoeff.h"
**
*     Input Variables
*
      double precision x,Q,y
**
*     Internal Variables
*
      double precision F2light,F2charm,F2bottom,F2top,F2total
      double precision FLlight,FLcharm,FLbottom,FLtop,FLtotal
      double precision F3light,F3charm,F3bottom,F3top,F3total
      double precision GetWMass,GetGFermi,GetProtonMass
      double precision GetSIATotalCrossSection
      double precision yp,ym,y2,ypc
      double precision Q2,MW,MW2,GF,GF2,MN
      double precision norm
      double precision conv
      character*21 obs
      double precision GetZMass,GetSin2ThetaW,xPDFj,AlphaQED ! EFT
      double precision e2,MZ2,xu,xubar,xd,xdbar,xs,xsbar,xc,xcbar
      double precision K,s2w,c2w  ! EFT
      double precision cF2u,cF2d,cF2u_d,cF2d_d,cF2s,cF2c ! EFT
      double precision cF3u,cF3d,cF3u_d,cF3d_d,cF3s,cF3c ! EFT
      parameter(conv=3.893793d10) ! conversion factor from GeV^-2 to 10^-38 cm^2
**
*     Output Variables
*
      double precision FKObservables
*
      obs = FKObservable
*
*     Useful definitions
*
      MW = GetWMass()
      MN = GetProtonMass()
      GF = GetGFermi()
*
      Q2  = Q * Q
      MW2 = MW * MW
      GF2 = GF * GF
      yp  = 1d0 + ( 1d0 - y )**2
      ym  = 1d0 - ( 1d0 - y )**2
      y2  = y * y
*
****  EFT cofficients
*
      e2 = 4d0 * pi * AlphaQED(Q)
c      au = 1.9d0/1d3**2
c      ad = 1.8d0/1d3**2
c      as =-3.0d0/1d3**2
c      ac = 5.0d0/1d3**2
cc     theory 181
c      au = 0.22d0/1d3**2
c      ad = -0.26d0/1d3**2
c      as = -0.26d0/1d3**2
c      ac = 0.22d0/1d3**2
cc     theory 182
c      au = -0.08d0/1d3**2
c      ad = +0.06d0/1d3**2
c      as = +0.06d0/1d3**2
c      ac = -0.08d0/1d3**2
cc     theory 183
c      au = +0.16d0/1d3**2
c      ad = 0.d0
c      as = 0.d0
c      ac = +0.16d0/1d3**2
cc     theory 184
c      au = -0.04d0/1d3**2
c      ad = -0.19d0/1d3**2
c      as = -0.19d0/1d3**2
c      ac = -0.04d0/1d3**2
c     theory 185
c      au = +0.57d0/1d3**2
c      ad = -0.48d0/1d3**2
c      as = +0.48d0/1d3**2
c      ac = -0.57d0/1d3**2
c     theory 186
c      au = -0.61d0/1d3**2
c      ad = +0.43d0/1d3**2
c      as = -0.43d0/1d3**2
c      ac = +0.61d0/1d3**2
c     theory 187
c      au = +0.15d0/1d3**2
c      ad = +0.09d0/1d3**2
c      as = -0.09d0/1d3**2
c      ac = -0.15d0/1d3**2
c     theory 188
c      au = 0.d0
c      ad = +0.11d0/1d3**2
c      as = -0.11d0/1d3**2
c      ac = 0.d0
c     theory 189
c      au = 0.33d0/1d3**2
c      ad = -0.33d0/1d3**2
c      as = 0.33d0/1d3**2
c      ac = -0.33d0/1d3**2
c     theory 190
c      au = -0.495d0/1d3**2
c      ad = 0.495d0/1d3**2
c      as = -0.495d0/1d3**2
c      ac = 0.495d0/1d3**2
c
      MZ2 = GetZMass()**2
      s2w = GetSin2ThetaW()
      c2w = 1d0 - s2w
      K = Q2 / 4d0 / c2w / s2w / (Q2+MZ2)
      xd   = xPDFj(1,x)
      xdbar= xPDFj(-1,x)
      xu   = xPDFj(2,x)
      xubar= xPDFj(-2,x)
      xs   = xPDFj(3,x)
      xsbar= xPDFj(-3,x)
      xc   = xPDFj(4,x)
      xcbar= xPDFj(-4,x)
      
******************************************************
*     F2 structure functions - extra EFT contributions
******************************************************
*
*     Up quark contribution
*      
      cF2u = 1d0 / 12d0 / e2**2  
     1     * (3d0*au*au*Q2*Q2 + 4d0*au*e2*Q2*(1d0+4d0*K*s2w*s2w))
     2     * (xu+xubar)
*     For deuteron data swap u <--> d
      cF2u_d = 1d0 / 24d0 / e2**2  
     1     * (3d0*au*au*Q2*Q2 + 4d0*au*e2*Q2*(1d0+4d0*K*s2w*s2w))
     2     * (xd+xdbar + xu+xubar)
*
*     Down quark contribution
*
      cF2d = 1d0 / 12d0 / e2**2 
     1     * (3d0*ad*ad*Q2*Q2 - 2d0*ad*e2*Q2*(1d0+4d0*K*s2w*s2w))
     2     * (xd+xdbar)
*     For deuteron data swap u <--> d
      cF2d_d = 1d0 / 24d0 / e2**2 
     1     * (3d0*ad*ad*Q2*Q2 - 2d0*ad*e2*Q2*(1d0+4d0*K*s2w*s2w))
     2     * (xu+xubar + xd+xdbar)
*
*     Strange
*
      cF2s = 1d0 / 12d0 / e2**2  
     1     * (3d0*as*as*Q2*Q2 - 2d0*as*e2*Q2*(1d0+4d0*K*s2w*s2w))
     2     * (xs+xsbar)
*
*     Charm
*
      cF2c = 1d0 / 12d0 / e2**2 
     1     * (3d0*ac*ac*Q2*Q2 + 4d0*ac*e2*Q2*(1d0+4d0*K*s2w*s2w))
     2     * (xc+xcbar)

******************************************************
*     xF3 structure function - extra EFT contributions
*     Note that cF3 give x*F3 to match the definition of 
*     F3light(x) in APFEL
******************************************************
*
*     Up quark contribution
*      
      cF3u = -1d0 / 12d0 / e2**2 !/x
     1     * (3d0*au*au*Q2*Q2 + 4d0*au*e2*Q2*(1d0+4d0*K*s2w*s2w))
     2     * (xu-xubar)
*     For deuteron data swap u <--> d
      cF3u_d = -1d0 / 24d0 / e2**2 !/x
     1     * (3d0*au*au*Q2*Q2 + 4d0*au*e2*Q2*(1d0+4d0*K*s2w*s2w))
     2     * (xu -xubar + xd-xdbar)
*
*     Down quark contribution
*      
      cF3d =  -1d0 / 12d0 / e2**2 !/x
     1     * (3d0*ad*ad*Q2*Q2 - 2d0*ad*e2*Q2*(1d0+4d0*K*s2w*s2w))
     2     * (xd-xdbar)
*     For deuteron data swap u <--> d
      cF3d_d =  -1d0 / 24d0 / e2**2 !/x
     1     * (3d0*ad*ad*Q2*Q2 - 2d0*ad*e2*Q2*(1d0+4d0*K*s2w*s2w))
     2     * (xu-xubar + xd-xdbar)
*
*     Strange
*
      cF3s =  -1d0 / 12d0 / e2**2!/x
     1     * (3d0*as*as*Q2*Q2 - 2d0*as*e2*Q2*(1d0+4d0*K*s2w*s2w))
     2     * (xs-xsbar)
*
*     Charm
*
      cF3c = -1d0 / 12d0 / e2**2!/x 
     1     *  (3d0*ac*ac*Q2*Q2 + 4d0*ac*e2*Q2*(1d0+4d0*K*s2w*s2w))
     2     *  (xc-xcbar)

******************************************************
*     Add them to explicit definition of observables
******************************************************
*     
****  Light structure function F2light
*
      if(obs(1:7).eq."DIS_F2L")then         
         FKObservables = F2light(x) + cF2u + cF2d + cF2s
*
****  Up structure function F2u
*
      elseif(obs(1:7).eq."DIS_F2U")then
         FKObservables = F2light(x) + cF2u
*
****  Down structure function F2d
*
      elseif(obs(1:7).eq."DIS_F2d")then
         FKObservables = F2light(x) + cF2d
*
****  Strange structure function F2s
*
      elseif(obs(1:7).eq."DIS_F2S")then
         FKObservables = F2light(x) + cF2s
*
****  Charm structure function F2charm
*
      elseif(obs(1:7).eq."DIS_F2C")then
         FKObservables = F2charm(x) + cF2c
*
****  Bottom structure function F2bottom
*
      elseif(obs(1:7).eq."DIS_F2B")then
         FKObservables = F2bottom(x)  ! no EFT contribution
*
****  Top structure function F2top
*
      elseif(obs(1:7).eq."DIS_F2T")then
         FKObservables = F2top(x)     ! no EFT contribution
*
****  Proton structure function F2 (neutral current)
*
      elseif(obs(1:10).eq."DIS_F2P_NC")then
         FKObservables = F2total(x) + cF2u + cF2d + cF2s + cF2c
*
****  Proton structure function F2 (electromagnetic) -- Do we use it in FK generation? 
*
      elseif(obs(1:7).eq."DIS_F2P")then
         FKObservables = F2total(x) ! no EFT contribution
*
****  Deuteron structure function F2  
*
*
      elseif(obs(1:7).eq."DIS_F2D")then
         FKObservables = F2total(x) + cF2u_d + cF2d_d + cF2s + cF2c
*
****  Light structure function FLlight
*
      elseif(obs(1:7).eq."DIS_FLL")then
         FKObservables = FLlight(x) ! no EFT contribution
*
****  Charm structure function FLcharm
*
      elseif(obs(1:7).eq."DIS_FLC")then
         FKObservables = FLcharm(x) ! no EFT contribution
*
****  Bottom structure function FLbottom
*
      elseif(obs(1:7).eq."DIS_FLB")then
         FKObservables = FLbottom(x) ! no EFT contribution
*
****  Top structure function FLtop
*
      elseif(obs(1:7).eq."DIS_FLT")then
         FKObservables = FLtop(x) ! no EFT contribution
*
****  Deuteron structure function FL
*
      elseif(obs(1:7).eq."DIS_FLD")then
         FKObservables = FLtotal(x)
*
****  Proton structure function FL (Neutral current)
*
      elseif(obs(1:10).eq."DIS_FLP_NC".or.
     1       obs(1:14).eq."DIS_FLP_CON_NC")then
         FKObservables = FLtotal(x)  ! no EFT contribution
*
****  Proton structure function FL (electromagnetic)
*
      elseif(obs(1:7).eq."DIS_FLP")then
         FKObservables = FLtotal(x)  ! no EFT contribution
*
****  F3 structure function
*
      elseif(obs(1:10).eq."DIS_F3P_NC")then
         FKObservables = F3total(x) + cF3u + cF3d + cF3s + cF3c
*
****  Electron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_NCE_L")then
         FKObservables = ( F2light(x) + cF2u + cF2d + cF2s
     1                 - ( y2 / yp ) * FLlight(x)
     2                 + ( ym / yp ) *(F3light(x) + cF3u + cF3d + cF3s))
*
****  Positron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_NCP_L")then
         FKObservables = ( F2light(x) + cF2u + cF2d + cF2s
     1                 - ( y2 / yp ) * FLlight(x)
     2                 - ( ym / yp ) *(F3light(x) + cF3u + cF3d + cF3s))
*
****  Electron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:10).eq."DIS_NCE_CH")then
         FKObservables = ( F2charm(x) + cF2c
     1                 - ( y2 / yp ) * FLcharm(x)
     2                 + ( ym / yp ) * (F3charm(x) + cF3c) )
*
****  Positron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:10).eq."DIS_NCP_CH")then
         FKObservables = ( F2charm(x) + cF2c
     1                 - ( y2 / yp ) * FLcharm(x)
     2                 - ( ym / yp ) * (F3charm(x) + cF3c) )
*
****  Electron scattering Reduced Cross-Section (bottom)
*
      elseif(obs(1:10).eq."DIS_NCE_BT")then  ! no EFT contribution
         FKObservables = ( F2bottom(x) 
     1                 - ( y2 / yp ) * FLbottom(x)
     2                 + ( ym / yp ) * F3bottom(x))
*
****  Positron scattering Reduced Cross-Section (bottom)
*
      elseif(obs(1:10).eq."DIS_NCP_BT")then  ! no EFT contribution
         FKObservables = ( F2bottom(x) 
     1                 - ( y2 / yp ) * FLbottom(x)
     2                 - ( ym / yp ) * F3bottom(x) )
*
****  Electron scattering Reduced Cross-Section (top)
*
      elseif(obs(1:10).eq."DIS_NCE_TP")then ! no EFT contribution
         FKObservables = ( F2top(x) 
     1                 - ( y2 / yp ) * FLtop(x)
     2                 + ( ym / yp ) * F3top(x) )
*
****  Positron scattering Reduced Cross-Section (top)
*
      elseif(obs(1:10).eq."DIS_NCP_TP")then ! no EFT contribution
         FKObservables = ( F2top(x) 
     1                 - ( y2 / yp ) * FLtop(x)
     2                 - ( ym / yp ) * F3top(x) )
*
****  Electron scattering Reduced Cross-Section on deuteron (inclusive)
*
      elseif(obs(1:9).eq."DIS_NCE_D")then
         FKObservables = ( F2total(x) + cF2u_d + cF2d_d + cF2s + cF2c
     1                 - ( y2 / yp ) * FLtotal(x)
     2        + ( ym / yp )
     3        * (F3total(x) + cF3u_d + cF3d_d + cF3s + cF3c) )
*
****  Positron scattering Reduced Cross-Section on deuteron (inclusive)
*
      elseif(obs(1:9).eq."DIS_NCP_D")then
         FKObservables = ( F2total(x) + cF2u_d + cF2d_d + cF2s + cF2c
     1                 - ( y2 / yp ) * FLtotal(x)
     2        - ( ym / yp )
     3        * (F3total(x) + cF3u_d + cF3d_d + cF3s + cF3c) )
*
****  Electron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_NCE")then
         FKObservables = ( F2total(x) + cF2u + cF2d + cF2s + cF2c
     1                 - ( y2 / yp ) * FLtotal(x)
     2        + ( ym / yp ) 
     3        * (F3total(x) + cF3u + cF3d + cF3s + cF3c) )
*
****  Positron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_NCP")then
         FKObservables = ( F2total(x) + cF2u + cF2d + cF2s + cF2c
     1                 - ( y2 / yp ) * FLtotal(x)
     2        - ( ym / yp )
     3        * (F3total(x) + cF3u + cF3d + cF3s + cF3c) )
*
****  Electron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_CCE_L")then ! no EFT contribution
         FKObservables = ( yp * F2light(x)
     1                 -   y2 * FLlight(x)
     2                 +   ym * F3light(x) )
         FKObservables = FKObservables / 4d0
*
****  Positron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_CCP_L")then ! no EFT contribution
         FKObservables = ( yp * F2light(x)
     1                 -   y2 * FLlight(x)
     2                 -   ym * F3light(x) )
         FKObservables = FKObservables / 4d0
*
****  Electron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_CCE_C")then ! no EFT contribution
         FKObservables = ( yp * F2charm(x)
     1                 -   y2 * FLcharm(x)
     2                 +   ym * F3charm(x) )
         FKObservables = FKObservables / 4d0
*
****  Positron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_CCP_C")then ! no EFT contribution
         FKObservables = ( yp * F2charm(x)  
     1                 -   y2 * FLcharm(x)
     2                 -   ym * F3charm(x)  )
         FKObservables = FKObservables / 4d0
*
****  Electron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_CCE")then! no EFT contribution
         FKObservables = ( yp * F2total(x)
     1                 -   y2 * FLtotal(x)
     2                 +   ym * F3total(x) )
         FKObservables = FKObservables / 4d0
*
****  Positron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_CCP")then ! no EFT contribution
         FKObservables = ( yp * F2total(x)
     1                 -   y2 * FLtotal(x)
     2                 -   ym * F3total(x) )
         FKObservables = FKObservables / 4d0
*
****  Neutrino scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_SNU_L")then ! no EFT contribution
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2light(x) 
     1                 -   y2  * FLlight(x)
     2                 +   ym  * F3light(x) )
         FKObservables = norm * FKObservables
*
****  Antineutrino scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_SNB_L")then ! no EFT contribution
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2light(x) 
     1                 -   y2  * FLlight(x)
     2                 -   ym  * F3light(x) )
         FKObservables = norm * FKObservables
*
****  Neutrino scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_SNU_C")then ! no EFT contribution
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2charm(x)
     1                 -   y2  * FLcharm(x)
     2                 +   ym  * F3charm(x) )
         FKObservables = norm * FKObservables
*
****  Antineutrino scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_SNB_C")then ! no EFT contribution
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2charm(x) 
     1                 -   y2  * FLcharm(x)
     2                 -   ym  * F3charm(x)  )
         FKObservables = norm * FKObservables
*
****  Neutrino scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_SNU")then ! no EFT contribution
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2total(x)
     1                 -   y2  * FLtotal(x)
     2                 +   ym  * F3total(x) )
         FKObservables = norm * FKObservables
*
****  Antineutrino scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_SNB")then ! no EFT contribution
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2total(x)
     1                 -   y2  * FLtotal(x)
     2                 -   ym  * F3total(x) )
         FKObservables = norm * FKObservables
*
****  Neutrino scattering Reduced Cross-Section (light)
*
      elseif(obs(1:12).eq."DIS_SNU_L_Pb")then ! no EFT contribution
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2light(x) 
     1                 -   y2  * FLlight(x)
     2                 +   ym  * F3light(x) )
         FKObservables = norm * FKObservables
*
****  Antineutrino scattering Reduced Cross-Section (light)
*
      elseif(obs(1:12).eq."DIS_SNB_L_Pb")then ! no EFT contribution
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2light(x) 
     1                 -   y2  * FLlight(x)
     2                 -   ym  * F3light(x) )
         FKObservables = norm * FKObservables
*
****  Neutrino scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:12).eq."DIS_SNU_C_Pb")then ! no EFT contribution
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2charm(x) 
     1                 -   y2  * FLcharm(x)
     2                 +   ym  * F3charm(x) )
         FKObservables = norm * FKObservables
*
****  Antineutrino scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:12).eq."DIS_SNB_C_Pb")then ! no EFT contribution
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2charm(x) 
     1                 -   y2  * FLcharm(x)
     2                 -   ym  * F3charm(x) )
         FKObservables = norm * FKObservables
*
****  Neutrino scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:10).eq."DIS_SNU_Pb")then ! no EFT contribution
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2total(x) 
     1                 -   y2  * FLtotal(x)
     2                 +   ym  * F3total(x) )
         FKObservables = norm * FKObservables
*
****  Antineutrino scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:10).eq."DIS_SNB_Pb")then ! no EFT contribution
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2total(x) 
     1                 -   y2  * FLtotal(x)
     2                 -   ym  * F3total(x) )
         FKObservables = norm * FKObservables
*
****  Dimuon neutrino and anti-neutrino cross section
*
************************************************************************
*
*     Note that NuTev dimuon data is given by the generator as:
*     
*                PI          dsigma^{2\mu}             1     1
*     100 * ------------ *  --------------- (E,x,y) * --- * ---   
*           G^2 M_p E_nu         dx dy                BrC   Acc
*     
*     Therefore, experimental acceptances and branching ratios are
*     corrected at the level of experimental data.
*     So we can compare directly with charm production expressions:
*
*         dsigma^{2\mu}            1     1        dsigma^{\nu, c}
*        -------------- (E,x,y) * --- * --- =   ----------------- 
*           dx dy                 BrC   Acc           dx dy
*
*     So by comparing with Eq. (116) of the "Notes on physical observables"
*     we find that one needs to multiply the output of the SIGMA_DM below 
*     by an overall factor:
*
*        dsigma^{\nu, c} |                 100        dsigma^{\nu, c} |
*        --------------- |        = --------------- * --------------- | 
*             dx dy      |_{data}   2*(1+Q2/MW2)**2        dx dy      |_{code}
*     
************************************************************************
*
****  Dimuon neutrino cross section
*
      elseif(obs(1:9).eq."DIS_DM_NU".or.obs(1:11).eq."DIS_DMN_CON")then ! no EFT contribution
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = 100d0 / 2d0 / ( 1d0 + Q2 / MW2 )**2 
         FKObservables = ( ypc * F2charm(x) 
     1                 -   y2  * FLcharm(x)
     2                 +   ym  * F3charm(x)  )
         FKObservables = norm * FKObservables
*
****  Dimuon anti-neutrino cross section
*
      elseif(obs(1:9).eq."DIS_DM_NB")then ! no EFT contribution
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = 100d0 / 2d0 / ( 1d0 + Q2 / MW2 )**2 
         FKObservables = ( ypc * F2charm(x) 
     1                 -   y2  * FLcharm(x)
     2                 -   ym  * F3charm(x) )
         FKObservables = norm * FKObservables
*
*     Single-Inclusive electron-positron annihilation (SIA)
*
****  SIA structure function F2 =  FT + FL
*
      elseif(obs(1:6).eq."SIA_F2")then
         FKObservables = F2total(x)
*
****  SIA structure function FL
*
      elseif(obs(1:6).eq."SIA_FL")then
         FKObservables = FLtotal(x)
*
****  SIA structure function FA
*
      elseif(obs(1:6).eq."SIA_FA")then
         FKObservables = F3total(x)
*
****  SIA absolute cross section (nf=4)
*
      elseif(obs(1:12).eq."SIA_XSEC_NF4")then
         FKObservables = GetSIATotalCrossSection(0,Q,"total")
     2                 * ( F2light(x) + F2charm(x) )
*
****  SIA absolute cross section
*
      elseif(obs(1:8).eq."SIA_XSEC")then
         FKObservables = GetSIATotalCrossSection(0,Q,"total")
     1                 * F2total(x)
*
****  SIA normalized light longitudinal cross section
*
      elseif(obs(1:20).eq."SIA_NORM_XSEC_LONG_L")then
         FKObservables = GetSIATotalCrossSection(0,Q,"total")
     1                 * FLlight(x)
     2                 / GetSIATotalCrossSection(ipt,Q,"light")
*
****  SIA normalized bottom longitudinal cross section
*
      elseif(obs(1:21).eq."SIA_NORM_XSEC_LONG_BT")then
         FKObservables = GetSIATotalCrossSection(0,Q,"total")
     1                 * FLbottom(x)
     2                 / GetSIATotalCrossSection(ipt,Q,"bottom")
*
****  SIA normalized total longitudinal cross section
*
      elseif(obs(1:18).eq."SIA_NORM_XSEC_LONG")then
         FKObservables = GetSIATotalCrossSection(0,Q,"total")
     1                 * FLtotal(x)
     2                 / GetSIATotalCrossSection(ipt,Q,"total")
*
****  SIA normalized light cross section
*
      elseif(obs(1:15).eq."SIA_NORM_XSEC_L")then
         FKObservables = GetSIATotalCrossSection(0,Q,"total")
     1                 * F2light(x)
     2                 / GetSIATotalCrossSection(ipt,Q,"light")
*
****  SIA normalized charm cross section
*
      elseif(obs(1:16).eq."SIA_NORM_XSEC_CH")then
         FKObservables = GetSIATotalCrossSection(0,Q,"total")
     1                 * F2charm(x)
     2                 / GetSIATotalCrossSection(ipt,Q,"charm")
*
****  SIA normalized bottom cross section
*
      elseif(obs(1:16).eq."SIA_NORM_XSEC_BT")then
         FKObservables = GetSIATotalCrossSection(0,Q,"total")
     1                 * F2bottom(x)
     2                 / GetSIATotalCrossSection(ipt,Q,"bottom")
*
****  SIA normalized top cross section
*
      elseif(obs(1:16).eq."SIA_NORM_XSEC_TP")then
         FKObservables = GetSIATotalCrossSection(0,Q,"total")
     1                 * F2top(x)
     2                 / GetSIATotalCrossSection(ipt,Q,"top")
*

****  SIA normalized cross section (nf=4)
*
      elseif(obs(1:17).eq."SIA_NORM_XSEC_NF4")then
         FKObservables = GetSIATotalCrossSection(0,Q,"total")
     1                 * ( F2light(x) + F2charm(x) )
     2                 / ( GetSIATotalCrossSection(ipt,Q,"light")
     3                   + GetSIATotalCrossSection(ipt,Q,"charm") )
*
****  SIA normalized total cross section
*
      elseif(obs(1:13).eq."SIA_NORM_XSEC")then
         FKObservables = GetSIATotalCrossSection(0,Q,"total")
     1                 * F2total(x)
     2                 / GetSIATotalCrossSection(ipt,Q,"total")
      else
         write(6,*) "In FKObservables.f:"
         write(6,*) "Invalid observable, obs = ",obs
         write(6,*) "  "
         call exit(-10)
      endif
*
      return
      end
