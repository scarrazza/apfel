************************************************************************
*
*     SetFKObservable.f:
*
*     This subroutine sets the observable according to the FKgenerator
*     naming.
*
************************************************************************
      subroutine SetFKObservable(obs)
*
      implicit none
*
      include "../commons/FKObservable.h"
**
*     Input Variables
*
      character*21 obs
*
      FKObservable = obs
*
****  PDFs for positivity observables
*
      if(obs(1:7).eq."DIS_XUQ".or.
     1     obs(1:7).eq."DIS_XUB".or.
     2     obs(1:7).eq."DIS_XDQ".or.
     3     obs(1:7).eq."DIS_XDB".or.
     4     obs(1:7).eq."DIS_XSQ".or.
     5     obs(1:7).eq."DIS_XSB".or.
     6     obs(1:7).eq."DIS_XCQ".or.
     7     obs(1:7).eq."DIS_XGL".or.
     8     obs(1:7).eq."DIS_XT3".or.
     9     obs(1:7).eq."DIS_XT8".or.
     1     obs(1:8).eq."DIS_XT15".or.
     2     obs(1:7).eq."DIS_XV3".or.
     3     obs(1:7).eq."DIS_XV8".or.
     4     obs(1:6).eq."DIS_XV")then
         call SetTargetDIS("proton") 
*
****  Light structure function F2light
*
      elseif(obs(1:7).eq."DIS_F2L")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Up structure function F2u
*
      elseif(obs(1:7).eq."DIS_F2U")then
         call SelectCharge("up")
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Down structure function F2d
*
      elseif(obs(1:7).eq."DIS_F2d")then
         call SelectCharge("down")
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Strange structure function F2s
*
      elseif(obs(1:7).eq."DIS_F2S")then
         call SelectCharge("strange")
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Proton structure function F2 (charm, CC, W-)
*
      elseif(obs(1:11).eq."DIS_F2C_CCE")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Proton structure function F2 (charm, CC, W+)
*
      elseif(obs(1:11).eq."DIS_F2C_CCP")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("proton")  
*
****  Charm structure function F2charm
*
      elseif(obs(1:7).eq."DIS_F2C")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")        
*
****  Bottom structure function F2bottom
*
      elseif(obs(1:7).eq."DIS_F2B")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Top structure function F2top
*
      elseif(obs(1:7).eq."DIS_F2T")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Deuteron structure function F2
*
      elseif(obs(1:7).eq."DIS_F2D")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("isoscalar")
*
****  Light structure function FLlight
*
      elseif(obs(1:7).eq."DIS_FLL")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Charm structure function FLcharm
*
      elseif(obs(1:7).eq."DIS_FLC")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Bottom structure function FLbottom
*
      elseif(obs(1:7).eq."DIS_FLB")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Top structure function FLtop
*
      elseif(obs(1:7).eq."DIS_FLT")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Deuteron structure function FL
*
      elseif(obs(1:7).eq."DIS_FLD")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("isoscalar")
*
****  Proton structure function F2 (Neutral current)
*
      elseif(obs(1:10).eq."DIS_F2P_NC")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Proton structure function F2 (electromagnetic)
*
      elseif(obs(1:7).eq."DIS_F2P")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")       
*
****  Proton structure function FL (Neutral current)
*
      elseif(obs(1:10).eq."DIS_FLP_NC".or.
     1       obs(1:14).eq."DIS_FLP_CON_NC")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Proton structure function FL (electromagnetic)
*
      elseif(obs(1:7).eq."DIS_FLP")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  F3 structure function
*
      elseif(obs(1:10).eq."DIS_F3P_NC")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")         
*
****  Electron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_NCE_L")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Positron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_NCP_L")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("proton")
*
****  Electron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:10).eq."DIS_NCE_CH")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Positron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:10).eq."DIS_NCP_CH")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("proton")
*
****  Electron scattering Reduced Cross-Section (bottom)
*
      elseif(obs(1:10).eq."DIS_NCE_BT")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Positron scattering Reduced Cross-Section (bottom)
*
      elseif(obs(1:10).eq."DIS_NCP_BT")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("proton")
*
****  Electron scattering Reduced Cross-Section (top)
*
      elseif(obs(1:10).eq."DIS_NCE_TP")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Positron scattering Reduced Cross-Section (top)
*
      elseif(obs(1:10).eq."DIS_NCP_TP")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("proton")
*
****  Electron scattering Reduced Cross-Section on deuteron (inclusive)
*
      elseif(obs(1:9).eq."DIS_NCE_D")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("isoscalar")
*
****  Positron scattering Reduced Cross-Section on deuteron (inclusive)
*
      elseif(obs(1:9).eq."DIS_NCP_D")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("isoscalar")
*
****  Electron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_NCE")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Positron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_NCP")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("proton")
*
****  Electron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_CCE_L")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Positron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_CCP_L")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("proton")
*
****  Electron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_CCE_C")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Positron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_CCP_C")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("proton")
*
****  Electron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_CCE")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Positron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_CCP")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("proton")
*
****  Neutrino scattering Reduced Cross-Section (light)
*
      elseif(obs(1:12).eq."DIS_SNU_L_Pb")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("neutrino")
         call SetTargetDIS("lead")
*
****  Antineutrino scattering Reduced Cross-Section (light)
*
      elseif(obs(1:12).eq."DIS_SNB_L_Pb")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("antineutrino")
         call SetTargetDIS("lead")
*
****  Neutrino scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:12).eq."DIS_SNU_C_Pb")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("neutrino")
         call SetTargetDIS("lead")
*
****  Antineutrino scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:12).eq."DIS_SNB_C_Pb")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("antineutrino")
         call SetTargetDIS("lead")
*
****  Neutrino scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:10).eq."DIS_SNU_Pb")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("neutrino")
         call SetTargetDIS("lead")
*
****  Antineutrino scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:10).eq."DIS_SNB_Pb")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("antineutrino")
         call SetTargetDIS("lead")
*
****  Neutrino scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_SNU_L")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("neutrino")
         call SetTargetDIS("isoscalar")
*
****  Antineutrino scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_SNB_L")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("antineutrino")
         call SetTargetDIS("isoscalar")
*
****  Neutrino scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_SNU_C")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("neutrino")
         call SetTargetDIS("isoscalar")
*
****  Antineutrino scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_SNB_C")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("antineutrino")
         call SetTargetDIS("isoscalar")
*
****  Neutrino scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_SNU")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("neutrino")
         call SetTargetDIS("isoscalar")
*
****  Antineutrino scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_SNB")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("antineutrino")
         call SetTargetDIS("isoscalar")
*
****  Dimuon neutrino cross section
*
      elseif(obs(1:9).eq."DIS_DM_NU".or.obs(1:11).eq."DIS_DMN_CON")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("neutrino")
         call SetTargetDIS("iron")
*
****  Dimuon anti-neutrino cross section
*
      elseif(obs(1:9).eq."DIS_DM_NB")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("antineutrino")
         call SetTargetDIS("iron")
*
*     Single-Inclusive electron-positron annihilation (SIA)
*
****  SIA structure function F2 =  FT + FL
*
      elseif(obs(1:6).eq."SIA_F2")then
         call SetTimeLikeEvolution(.true.)
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  SIA structure function FL
*
      elseif(obs(1:6).eq."SIA_FL")then
         call SetTimeLikeEvolution(.true.)
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  SIA structure function FA
*
      elseif(obs(1:6).eq."SIA_FA")then
         call SetTimeLikeEvolution(.true.)
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  SIA absolute cross section (nf=4)
*
      elseif(obs(1:12).eq."SIA_XSEC_NF4")then
         call SetTimeLikeEvolution(.true.)
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  SIA absolute cross section
*
      elseif(obs(1:8).eq."SIA_XSEC")then
         call SetTimeLikeEvolution(.true.)
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  SIA normalized light longitudinal cross section
*
      elseif(obs(1:20).eq."SIA_NORM_XSEC_LONG_L")then
         call SetTimeLikeEvolution(.true.)
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  SIA normalized bottom longitudinal cross section
*
      elseif(obs(1:21).eq."SIA_NORM_XSEC_LONG_BT")then
         call SetTimeLikeEvolution(.true.)
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  SIA normalized total longitudinal cross section
*
      elseif(obs(1:18).eq."SIA_NORM_XSEC_LONG")then
         call SetTimeLikeEvolution(.true.)
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  SIA normalized light cross section
*
      elseif(obs(1:15).eq."SIA_NORM_XSEC_L")then
         call SetTimeLikeEvolution(.true.)
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  SIA normalized charm cross section
*
      elseif(obs(1:16).eq."SIA_NORM_XSEC_CH")then
         call SetTimeLikeEvolution(.true.)
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  SIA normalized bottom cross section
*
      elseif(obs(1:16).eq."SIA_NORM_XSEC_BT")then
         call SetTimeLikeEvolution(.true.)
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  SIA normalized top cross section
*
      elseif(obs(1:16).eq."SIA_NORM_XSEC_TP")then
         call SetTimeLikeEvolution(.true.)
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  SIA normalized total cross section (nf=4)
*
      elseif(obs(1:17).eq."SIA_NORM_XSEC_NF4")then
         call SetTimeLikeEvolution(.true.)
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  SIA normalized total cross section
*
      elseif(obs(1:13).eq."SIA_NORM_XSEC")then
         call SetTimeLikeEvolution(.true.)
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
      else
         write(6,*) "In SetFKObservable.f:"
         write(6,*) "Invalid observable, obs = ",obs
         write(6,*) "  "
         call exit(-10)
      endif
*
      return
      end
