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
      double precision xPDF
      double precision GetWMass,GetGFermi,GetProtonMass
      double precision GetSIATotalCrossSection
      double precision yp,ym,y2,ypc
      double precision Q2,MW,MW2,GF,GF2,MN
      double precision norm
      double precision conv
      character*21 obs
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
****  Up quark PDF
*
      if(obs(1:7).eq."DIS_XUQ")then
         FKObservables = xPDF(2,x)
*
****  Up antiquark PDF
*
      elseif(obs(1:7).eq."DIS_XUB")then
         FKObservables = xPDF(-2,x)
*
****  Down quark PDF
*
      elseif(obs(1:7).eq."DIS_XDQ")then
         FKObservables = xPDF(1,x)
*
****  Down antiquark PDF
*
      elseif(obs(1:7).eq."DIS_XDB")then
         FKObservables = xPDF(-1,x)
*
****  Strange quark PDF
*
      elseif(obs(1:7).eq."DIS_XSQ")then
         FKObservables = xPDF(3,x)
*
****  Strange antiquark PDF
*
      elseif(obs(1:7).eq."DIS_XSB")then
         FKObservables = xPDF(-3,x)
*
****  Charm quark PDF
*
      elseif(obs(1:7).eq."DIS_XCQ")then
         FKObservables = xPDF(4,x)
*
****  Gluon PDF
*
      elseif(obs(1:7).eq."DIS_XGL")then
         FKObservables = xPDF(0,x)
*
****  Non-singlet T3
*
      elseif(obs(1:7).eq."DIS_XT3")then
         FKObservables = xPDF(+2,x) + xpdf(-2,x)
     1        - 1d0 * ( xPDF(+1,x) + xPDF(-1,x) )
*
****  Non-singlet T8
*
      elseif(obs(1:7).eq."DIS_XT8")then
         FKObservables = xPDF(+2,x) + xpdf(-2,x)
     1        + 1d0 * ( xPDF(+1,x) + xPDF(-1,x) )
     1        - 2d0 * ( xPDF(+3,x) + xPDF(-3,x) )
*
****  Non-singlet T15
*
      elseif(obs(1:8).eq."DIS_XT15")then
         FKObservables = xPDF(+2,x) + xpdf(-2,x)
     1        + 1d0 * ( xPDF(+1,x) + xPDF(-1,x) )
     1        + 1d0 * ( xPDF(+3,x) + xPDF(-3,x) )
     1        - 3d0 * ( xPDF(+4,x) + xPDF(-4,x) )
*
****  Valence V3
*
      elseif(obs(1:7).eq."DIS_XV3")then
         FKObservables = xPDF(+2,x) - xPDF(-2,x) 
     1        - 1d0 * ( xPDF(+1,x) - xPDF(-1,x) )
*
****  Valence V8
*
      elseif(obs(1:7).eq."DIS_XV8")then
         FKObservables = xPDF(+2,x) - xPDF(-2,x) 
     1        + 1d0 * ( xPDF(+1,x) - xPDF(-1,x) )
     2        - 2d0 * ( xPDF(+3,x) - xPDF(-3,x) )
*
*
****  Valence V
*
      elseif(obs(1:6).eq."DIS_XV")then
         FKObservables = xPDF(+2,x) - xPDF(-2,x) 
     1        + xPDF(+1,x) - xPDF(-1,x)
     1        + xPDF(+3,x) - xPDF(-3,x)
     1        + xPDF(+4,x) - xPDF(-4,x) 
     1        + xPDF(+5,x) - xPDF(-5,x)
     1        + xPDF(+6,x) - xPDF(-6,x) 
*
****  Light structure function F2light
*
      elseif(obs(1:7).eq."DIS_F2L")then         
         FKObservables = F2light(x)
*
****  Up structure function F2u
*
      elseif(obs(1:7).eq."DIS_F2U")then
         FKObservables = F2light(x)
*
****  Down structure function F2d
*
      elseif(obs(1:7).eq."DIS_F2d")then
         FKObservables = F2light(x)
*
****  Strange structure function F2s
*
      elseif(obs(1:7).eq."DIS_F2S")then
         FKObservables = F2light(x)
*
****  Proton structure function F2 (charm, CC, W-)
*
      elseif(obs(1:11).eq."DIS_F2C_CCE")then
         FKObservables = F2charm(x)
*
****  Proton structure function F2 (charm, CC, W+)
*
      elseif(obs(1:11).eq."DIS_F2C_CCP")then
         FKObservables = F2charm(x)        
*
****  Charm structure function F2charm
*
      elseif(obs(1:7).eq."DIS_F2C")then
         FKObservables = F2charm(x)    
*
****  Bottom structure function F2bottom
*
      elseif(obs(1:7).eq."DIS_F2B")then
         FKObservables = F2bottom(x)
*
****  Top structure function F2top
*
      elseif(obs(1:7).eq."DIS_F2T")then
         FKObservables = F2top(x)
*
****  Proton structure function F2 (neutral current)
*
      elseif(obs(1:10).eq."DIS_F2P_NC")then
         FKObservables = F2total(x)
*
****  Proton structure function F2 (electromagnetic)
*
      elseif(obs(1:7).eq."DIS_F2P")then
         FKObservables = F2total(x)
*
****  Deuteron structure function F2
*
      elseif(obs(1:7).eq."DIS_F2D")then
         FKObservables = F2total(x)         
*
****  Light structure function FLlight
*
      elseif(obs(1:7).eq."DIS_FLL")then
         FKObservables = FLlight(x)
*
****  Charm structure function FLcharm
*
      elseif(obs(1:7).eq."DIS_FLC")then
         FKObservables = FLcharm(x)
*
****  Bottom structure function FLbottom
*
      elseif(obs(1:7).eq."DIS_FLB")then
         FKObservables = FLbottom(x)
*
****  Top structure function FLtop
*
      elseif(obs(1:7).eq."DIS_FLT")then
         FKObservables = FLtop(x)
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
         FKObservables = FLtotal(x)
*
****  Proton structure function FL (electromagnetic)
*
      elseif(obs(1:7).eq."DIS_FLP")then
         FKObservables = FLtotal(x)
*
****  F3 structure function
*
      elseif(obs(1:10).eq."DIS_F3P_NC")then
         FKObservables = F3total(x)           
*
****  Electron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_NCE_L")then
         FKObservables = ( F2light(x) 
     1                 - ( y2 / yp ) * FLlight(x)
     2                 + ( ym / yp ) * F3light(x) )
*
****  Positron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_NCP_L")then
         FKObservables = ( F2light(x) 
     1                 - ( y2 / yp ) * FLlight(x)
     2                 - ( ym / yp ) * F3light(x) )
*
****  Electron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:10).eq."DIS_NCE_CH")then
         FKObservables = ( F2charm(x) 
     1                 - ( y2 / yp ) * FLcharm(x)
     2                 + ( ym / yp ) * F3charm(x) )
*
****  Positron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:10).eq."DIS_NCP_CH")then
         FKObservables = ( F2charm(x) 
     1                 - ( y2 / yp ) * FLcharm(x)
     2                 - ( ym / yp ) * F3charm(x) )
*
****  Electron scattering Reduced Cross-Section (bottom)
*
      elseif(obs(1:10).eq."DIS_NCE_BT")then
         FKObservables = ( F2bottom(x) 
     1                 - ( y2 / yp ) * FLbottom(x)
     2                 + ( ym / yp ) * F3bottom(x) )
*
****  Positron scattering Reduced Cross-Section (bottom)
*
      elseif(obs(1:10).eq."DIS_NCP_BT")then
         FKObservables = ( F2bottom(x) 
     1                 - ( y2 / yp ) * FLbottom(x)
     2                 - ( ym / yp ) * F3bottom(x) )
*
****  Electron scattering Reduced Cross-Section (top)
*
      elseif(obs(1:10).eq."DIS_NCE_TP")then
         FKObservables = ( F2top(x) 
     1                 - ( y2 / yp ) * FLtop(x)
     2                 + ( ym / yp ) * F3top(x) )
*
****  Positron scattering Reduced Cross-Section (top)
*
      elseif(obs(1:10).eq."DIS_NCP_TP")then
         FKObservables = ( F2top(x) 
     1                 - ( y2 / yp ) * FLtop(x)
     2                 - ( ym / yp ) * F3top(x) )
*
****  Electron scattering Reduced Cross-Section on deuteron (inclusive)
*
      elseif(obs(1:9).eq."DIS_NCE_D")then
         FKObservables = ( F2total(x) 
     1                 - ( y2 / yp ) * FLtotal(x)
     2                 + ( ym / yp ) * F3total(x) )
*
****  Positron scattering Reduced Cross-Section on deuteron (inclusive)
*
      elseif(obs(1:9).eq."DIS_NCP_D")then
         FKObservables = ( F2total(x) 
     1                 - ( y2 / yp ) * FLtotal(x)
     2                 - ( ym / yp ) * F3total(x) )
*
****  Electron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_NCE")then
         FKObservables = ( F2total(x) 
     1                 - ( y2 / yp ) * FLtotal(x)
     2                 + ( ym / yp ) * F3total(x) )
*
****  Positron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_NCP")then
         FKObservables = ( F2total(x) 
     1                 - ( y2 / yp ) * FLtotal(x)
     2                 - ( ym / yp ) * F3total(x) )
*
****  Electron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_CCE_L")then
         FKObservables = ( yp * F2light(x) 
     1                 -   y2 * FLlight(x)
     2                 +   ym * F3light(x) )
         FKObservables = FKObservables / 4d0
*
****  Positron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_CCP_L")then
         FKObservables = ( yp * F2light(x) 
     1                 -   y2 * FLlight(x)
     2                 -   ym * F3light(x) )
         FKObservables = FKObservables / 4d0
*
****  Electron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_CCE_C")then
         FKObservables = ( yp * F2charm(x) 
     1                 -   y2 * FLcharm(x)
     2                 +   ym * F3charm(x) )
         FKObservables = FKObservables / 4d0
*
****  Positron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_CCP_C")then
         FKObservables = ( yp * F2charm(x) 
     1                 -   y2 * FLcharm(x)
     2                 -   ym * F3charm(x) )
         FKObservables = FKObservables / 4d0
*
****  Electron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_CCE")then
         FKObservables = ( yp * F2total(x) 
     1                 -   y2 * FLtotal(x)
     2                 +   ym * F3total(x) )
         FKObservables = FKObservables / 4d0
*
****  Positron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_CCP")then
         FKObservables = ( yp * F2total(x) 
     1                 -   y2 * FLtotal(x)
     2                 -   ym * F3total(x) )
         FKObservables = FKObservables / 4d0
*
****  Neutrino scattering Reduced Cross-Section (light)
*
      elseif(obs(1:12).eq."DIS_SNU_L_Pb")then
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2light(x) 
     1                 -   y2  * FLlight(x)
     2                 +   ym  * F3light(x) )
         FKObservables = norm * FKObservables
*
****  Antineutrino scattering Reduced Cross-Section (light)
*
      elseif(obs(1:12).eq."DIS_SNB_L_Pb")then
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2light(x) 
     1                 -   y2  * FLlight(x)
     2                 -   ym  * F3light(x) )
         FKObservables = norm * FKObservables
*
****  Neutrino scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:12).eq."DIS_SNU_C_Pb")then
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2charm(x) 
     1                 -   y2  * FLcharm(x)
     2                 +   ym  * F3charm(x) )
         FKObservables = norm * FKObservables
*
****  Antineutrino scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:12).eq."DIS_SNB_C_Pb")then
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2charm(x) 
     1                 -   y2  * FLcharm(x)
     2                 -   ym  * F3charm(x) )
         FKObservables = norm * FKObservables
*
****  Neutrino scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:10).eq."DIS_SNU_Pb")then
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2total(x) 
     1                 -   y2  * FLtotal(x)
     2                 +   ym  * F3total(x) )
         FKObservables = norm * FKObservables
*
****  Antineutrino scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:10).eq."DIS_SNB_Pb")then
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2total(x) 
     1                 -   y2  * FLtotal(x)
     2                 -   ym  * F3total(x) )
         FKObservables = norm * FKObservables
*
****  Neutrino scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_SNU_L")then
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2light(x) 
     1                 -   y2  * FLlight(x)
     2                 +   ym  * F3light(x) )
         FKObservables = norm * FKObservables
*
****  Antineutrino scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_SNB_L")then
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2light(x) 
     1                 -   y2  * FLlight(x)
     2                 -   ym  * F3light(x) )
         FKObservables = norm * FKObservables
*
****  Neutrino scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_SNU_C")then
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2charm(x) 
     1                 -   y2  * FLcharm(x)
     2                 +   ym  * F3charm(x) )
         FKObservables = norm * FKObservables
*
****  Antineutrino scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_SNB_C")then
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2charm(x) 
     1                 -   y2  * FLcharm(x)
     2                 -   ym  * F3charm(x) )
         FKObservables = norm * FKObservables
*
****  Neutrino scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_SNU")then
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = conv * GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2 )
         FKObservables = ( ypc * F2total(x) 
     1                 -   y2  * FLtotal(x)
     2                 +   ym  * F3total(x) )
         FKObservables = norm * FKObservables
*
****  Antineutrino scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_SNB")then
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
      elseif(obs(1:9).eq."DIS_DM_NU".or.obs(1:11).eq."DIS_DMN_CON")then
         ypc  = yp - 2d0 * ( MN * x * y )**2 / Q2
         norm = 100d0 / 2d0 / ( 1d0 + Q2 / MW2 )**2 
         FKObservables = ( ypc * F2charm(x) 
     1                 -   y2  * FLcharm(x)
     2                 +   ym  * F3charm(x) )
         FKObservables = norm * FKObservables
*
****  Dimuon anti-neutrino cross section
*
      elseif(obs(1:9).eq."DIS_DM_NB")then
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
