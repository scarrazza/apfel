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
**
*     Input Variables
*
      character*15 obs
**
*     Electromagnetic Observables
*
****  Light structure function F2light
*
      if(obs(1:7).eq."DIS_F2L")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Up structure function F2u
*
      elseif(obs(1:7).eq."DIS_F2U")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Down structure function F2d
*
      elseif(obs(1:7).eq."DIS_F2d")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  Strange structure function F2s
*
      elseif(obs(1:7).eq."DIS_F2S")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
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
****  Proton structure function F2
*
      elseif(obs(1:7).eq."DIS_F2P")then
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
****  Proton structure function FL
*
      elseif(obs(1:7).eq."DIS_FLP")then
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
**
*     Neutral-Current Observables
*
****  F2 structure function
*
      elseif(obs(1:10).eq."DIS_F2P_NC")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
*
****  FL longitudinal structure function
*
      elseif(obs(1:10).eq."DIS_FLP_NC".or.
     1       obs(1:14).eq."DIS_FLP_CON_NC")then
         call SetProcessDIS("NC")
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
**
*     Charged-Current Observables
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
**
*     Neutrino Observables
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
      else
         write(6,*) "In SetFKObservable.f:"
         write(6,*) "Invalid observable, obs = ",obs
         write(6,*) "  "
         call exit(-10)
      endif
*
      return
      end
