************************************************************************
*
*     FKSimulator.f:
*
*     This function returns simulates the FKgenerator of the NNPDF code.
*
************************************************************************
      function FKSimulator(obs,x,Q,y,i,beta)
*
      implicit none
*
      include "../commons/EvolOp.h"
      include "../commons/grid.h"
      include "../commons/DISOperators.h"
      include "../commons/WMass.h"
      include "../commons/ProtonMass.h"
      include "../commons/GFermi.h"
      include "../commons/consts.h"
**
*     Input Variables
*
      integer i,beta
      double precision x,Q,y
      character*15 obs
**
*     Internal Variables
*
      integer n
      integer alpha
      double precision w_int_gen
      double precision yp,ym,y2,ypc
      double precision Q2,MW2,GF2,MN
      double precision norm
      double precision conv
      parameter(conv=3.893793d10) ! conversion factor from GeV^-2 to 10^-38 cm^2
**
*     Output Variables
*
      double precision FKSimulator
*
*     Check whether the evolution operator has been actually computed
*
      if(.not.EvolOp)then
         write(6,*) "The evolution operator computation is disabled."
         write(6,*) "The 'FKSimulator' function cannot be used."
         write(6,*) "   "
         call exit(-10)
      endif
*
      if(i.lt.0.or.i.gt.13)then
         write(6,*) "In FKSimulator.f:"
         write(6,*) "Invalid index, i =",i
         call exit(-10)
      endif
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In FKSimulator.f:"
         write(6,*) "Invalid value of x =",x
         call exit(-10)
      endif
*
      if(beta.lt.0.or.beta.gt.nin(0))then
         write(6,*) "In FKSimulator.f:"
         write(6,*) "Invalid index, beta =",beta
         call exit(-10)
      endif
*
*     Useful definitions
*
      Q2  = Q * Q
      MW2 = MW * MW
      GF2 = GFermi * GFermi
      MN  = MProton
      yp  = 1d0 + ( 1d0 - y )**2d0
      ym  = 1d0 - ( 1d0 - y )**2d0
      y2  = y * y
*
*     Interpolate
*
      n = inter_degree(0)
      FKSimulator = 0d0
**
*     Electromagnetic Observables
*
****  Light structure function F2light
*
      if(obs(1:7).eq."DIS_F2L")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * EvOpF2(3,i,alpha,beta)
         enddo
*
****  Up structure function F2u
*
      elseif(obs(1:7).eq."DIS_F2U")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         write(6,*) "For this observables I don't know yet what to do!"
         call exit(-10)
*
****  Down structure function F2d
*
      elseif(obs(1:7).eq."DIS_F2d")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         write(6,*) "For this observables I don't know yet what to do!"
         call exit(-10)
*
****  Strange structure function F2s
*
      elseif(obs(1:7).eq."DIS_F2S")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         write(6,*) "For this observables I don't know yet what to do!"
         call exit(-10)
*
****  Charm structure function F2charm
*
      elseif(obs(1:7).eq."DIS_F2C")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * EvOpF2(4,i,alpha,beta)
         enddo
*
****  Bottom structure function F2bottom
*
      elseif(obs(1:7).eq."DIS_F2B")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * EvOpF2(5,i,alpha,beta)
         enddo
*
****  Top structure function F2top
*
      elseif(obs(1:7).eq."DIS_F2T")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * EvOpF2(6,i,alpha,beta)
         enddo
*
****  Proton structure function F2
*
      elseif(obs(1:7).eq."DIS_F2P")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * EvOpF2(7,i,alpha,beta)
         enddo
*
****  Deuteron structure function F2
*
      elseif(obs(1:7).eq."DIS_F2D")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("isoscalar")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * EvOpF2(7,i,alpha,beta)
         enddo
*
****  Light structure function FLlight
*
      elseif(obs(1:7).eq."DIS_FLL")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * EvOpFL(3,i,alpha,beta)
         enddo
*
****  Charm structure function FLcharm
*
      elseif(obs(1:7).eq."DIS_FLC")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * EvOpFL(4,i,alpha,beta)
         enddo
*
****  Bottom structure function FLbottom
*
      elseif(obs(1:7).eq."DIS_FLB")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * EvOpFL(5,i,alpha,beta)
         enddo
*
****  Top structure function FLtop
*
      elseif(obs(1:7).eq."DIS_FLT")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * EvOpFL(6,i,alpha,beta)
         enddo
*
****  Proton structure function FL
*
      elseif(obs(1:7).eq."DIS_FLP")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * EvOpFL(7,i,alpha,beta)
         enddo
*
****  Deuteron structure function FL
*
      elseif(obs(1:7).eq."DIS_FLD")then
         call SetProcessDIS("EM")
         call SetProjectileDIS("electron")
         call SetTargetDIS("isoscalar")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * EvOpFL(7,i,alpha,beta)
         enddo
**
*     Neutral-Current Observables
*
****  F2 structure function
*
      elseif(obs(1:10).eq."DIS_F2P_NC")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * EvOpF2(7,i,alpha,beta)
         enddo
*
****  FL longitudinal structure function
*
      elseif(obs(1:10).eq."DIS_FLP_NC".or.
     1       obs(1:14).eq."DIS_FLP_CON_NC")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * EvOpFL(7,i,alpha,beta)
         enddo
*
****  F3 structure function
*
      elseif(obs(1:10).eq."DIS_F3P_NC")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * EvOpF3(7,i,alpha,beta)
         enddo
*
****  Electron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_NCE")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( EvOpF2(7,i,alpha,beta) 
     2                  - ( y2 / yp ) * EvOpFL(7,i,alpha,beta)
     3                  + ( ym / yp ) * EvOpF3(7,i,alpha,beta) )
         enddo
*
****  Positron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_NCP")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( EvOpF2(7,i,alpha,beta) 
     2                  - ( y2 / yp ) * EvOpFL(7,i,alpha,beta)
     3                  - ( ym / yp ) * EvOpF3(7,i,alpha,beta) )
         enddo
*
****  Electron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_NCE_L")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( EvOpF2(3,i,alpha,beta) 
     2                  - ( y2 / yp ) * EvOpFL(3,i,alpha,beta)
     3                  + ( ym / yp ) * EvOpF3(3,i,alpha,beta) )
         enddo
*
****  Positron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_NCP_L")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( EvOpF2(3,i,alpha,beta) 
     2                  - ( y2 / yp ) * EvOpFL(3,i,alpha,beta)
     3                  - ( ym / yp ) * EvOpF3(3,i,alpha,beta) )
         enddo
*
****  Electron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:10).eq."DIS_NCE_CH")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( EvOpF2(4,i,alpha,beta) 
     2                  - ( y2 / yp ) * EvOpFL(4,i,alpha,beta)
     3                  + ( ym / yp ) * EvOpF3(4,i,alpha,beta) )
         enddo
*
****  Positron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:10).eq."DIS_NCP_CH")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( EvOpF2(4,i,alpha,beta) 
     2                  - ( y2 / yp ) * EvOpFL(4,i,alpha,beta)
     3                  - ( ym / yp ) * EvOpF3(4,i,alpha,beta) )
         enddo
*
****  Electron scattering Reduced Cross-Section (bottom)
*
      elseif(obs(1:10).eq."DIS_NCE_BT")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( EvOpF2(5,i,alpha,beta) 
     2                  - ( y2 / yp ) * EvOpFL(5,i,alpha,beta)
     3                  + ( ym / yp ) * EvOpF3(5,i,alpha,beta) )
         enddo
*
****  Positron scattering Reduced Cross-Section (bottom)
*
      elseif(obs(1:10).eq."DIS_NCP_BT")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( EvOpF2(5,i,alpha,beta) 
     2                  - ( y2 / yp ) * EvOpFL(5,i,alpha,beta)
     3                  - ( ym / yp ) * EvOpF3(5,i,alpha,beta) )
         enddo
*
****  Electron scattering Reduced Cross-Section (top)
*
      elseif(obs(1:10).eq."DIS_NCE_TP")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( EvOpF2(6,i,alpha,beta) 
     2                  - ( y2 / yp ) * EvOpFL(6,i,alpha,beta)
     3                  + ( ym / yp ) * EvOpF3(6,i,alpha,beta) )
         enddo
*
****  Positron scattering Reduced Cross-Section (top)
*
      elseif(obs(1:10).eq."DIS_NCP_TP")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( EvOpF2(6,i,alpha,beta) 
     2                  - ( y2 / yp ) * EvOpFL(6,i,alpha,beta)
     3                  - ( ym / yp ) * EvOpF3(6,i,alpha,beta) )
         enddo
*
****  Electron scattering Reduced Cross-Section on deuteron (inclusive)
*
      elseif(obs(1:9).eq."DIS_NCE_D")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("isoscalar")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( EvOpF2(7,i,alpha,beta) 
     2                  - ( y2 / yp ) * EvOpFL(7,i,alpha,beta)
     3                  + ( ym / yp ) * EvOpF3(7,i,alpha,beta) )
         enddo
*
****  Positron scattering Reduced Cross-Section on deuteron (inclusive)
*
      elseif(obs(1:9).eq."DIS_NCP_D")then
         call SetProcessDIS("NC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("isoscalar")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( EvOpF2(7,i,alpha,beta) 
     2                  - ( y2 / yp ) * EvOpFL(7,i,alpha,beta)
     3                  - ( ym / yp ) * EvOpF3(7,i,alpha,beta) )
         enddo
**
*     Charged-Current Observables
*
****  Electron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_CCE")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( yp * EvOpF2(7,i,alpha,beta) 
     2                  -   y2 * EvOpFL(7,i,alpha,beta)
     3                  -   ym * EvOpF3(7,i,alpha,beta) )
         enddo
         FKSimulator = FKSimulator / 4d0
*
****  Positron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_CCP")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( yp * EvOpF2(7,i,alpha,beta) 
     2                  -   y2 * EvOpFL(7,i,alpha,beta)
     3                  +   ym * EvOpF3(7,i,alpha,beta) )
         enddo
         FKSimulator = FKSimulator / 4d0
*
****  Electron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_CCE_L")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( yp * EvOpF2(3,i,alpha,beta) 
     2                  -   y2 * EvOpFL(3,i,alpha,beta)
     3                  -   ym * EvOpF3(3,i,alpha,beta) )
         enddo
         FKSimulator = FKSimulator / 4d0
*
****  Positron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_CCP_L")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( yp * EvOpF2(3,i,alpha,beta) 
     2                  -   y2 * EvOpFL(3,i,alpha,beta)
     3                  +   ym * EvOpF3(3,i,alpha,beta) )
         enddo
         FKSimulator = FKSimulator / 4d0
*
****  Electron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_CCE_C")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("electron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( yp * EvOpF2(4,i,alpha,beta) 
     2                  -   y2 * EvOpFL(4,i,alpha,beta)
     3                  -   ym * EvOpF3(4,i,alpha,beta) )
         enddo
         FKSimulator = FKSimulator / 4d0
*
****  Positron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_CCP_C")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("positron")
         call SetTargetDIS("proton")
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( yp * EvOpF2(4,i,alpha,beta) 
     2                  -   y2 * EvOpFL(4,i,alpha,beta)
     3                  +   ym * EvOpF3(4,i,alpha,beta) )
         enddo
         FKSimulator = FKSimulator / 4d0
**
*     Neutrino Observables
*
****  Neutrino scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_SNU")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("neutrino")
         call SetTargetDIS("isoscalar")
         ypc  = yp - 2d0 * ( MN * x * y )**2d0 / Q2
         norm = GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2d0 )
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( ypc * EvOpF2(7,i,alpha,beta) 
     2                  -   y2  * EvOpFL(7,i,alpha,beta)
     3                  +   ym  * EvOpF3(7,i,alpha,beta) )
         enddo
         FKSimulator = norm * FKSimulator
*
****  Antineutrino scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_SNB")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("antineutrino")
         call SetTargetDIS("isoscalar")
         ypc  = yp - 2d0 * ( MN * x * y )**2d0 / Q2
         norm = GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2d0 )
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( ypc * EvOpF2(7,i,alpha,beta) 
     2                  -   y2  * EvOpFL(7,i,alpha,beta)
     3                  -   ym  * EvOpF3(7,i,alpha,beta) )
         enddo
         FKSimulator = norm * FKSimulator
*
****  Neutrino scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_SNU_L")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("neutrino")
         call SetTargetDIS("isoscalar")
         ypc  = yp - 2d0 * ( MN * x * y )**2d0 / Q2
         norm = GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2d0 )
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( ypc * EvOpF2(3,i,alpha,beta) 
     2                  -   y2  * EvOpFL(3,i,alpha,beta)
     3                  +   ym  * EvOpF3(3,i,alpha,beta) )
         enddo
         FKSimulator = norm * FKSimulator
*
****  Antineutrino scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_SNB_L")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("antineutrino")
         call SetTargetDIS("isoscalar")
         ypc  = yp - 2d0 * ( MN * x * y )**2d0 / Q2
         norm = GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2d0 )
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( ypc * EvOpF2(3,i,alpha,beta) 
     2                  -   y2  * EvOpFL(3,i,alpha,beta)
     3                  -   ym  * EvOpF3(3,i,alpha,beta) )
         enddo
         FKSimulator = norm * FKSimulator
*
****  Neutrino scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_SNU_C")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("neutrino")
         call SetTargetDIS("isoscalar")
         ypc  = yp - 2d0 * ( MN * x * y )**2d0 / Q2
         norm = GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2d0 )
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( ypc * EvOpF2(4,i,alpha,beta) 
     2                  -   y2  * EvOpFL(4,i,alpha,beta)
     3                  +   ym  * EvOpF3(4,i,alpha,beta) )
         enddo
         FKSimulator = norm * FKSimulator
*
****  Antineutrino scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_SNB_C")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("antineutrino")
         call SetTargetDIS("isoscalar")
         ypc  = yp - 2d0 * ( MN * x * y )**2d0 / Q2
         norm = GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2d0 )
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( ypc * EvOpF2(4,i,alpha,beta) 
     2                  -   y2  * EvOpFL(4,i,alpha,beta)
     3                  -   ym  * EvOpF3(4,i,alpha,beta) )
         enddo
         FKSimulator = norm * FKSimulator
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
         call SetProcessDIS("CC")
         call SetProjectileDIS("neutrino")
         call SetTargetDIS("iron")
         ypc  = yp - 2d0 * ( MN * x * y )**2d0 / Q2
         norm = 100d0 / 2d0 / ( 1d0 + Q2 / MW2 )**2d0 
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( ypc * EvOpF2(4,i,alpha,beta) 
     2                  -   y2  * EvOpFL(4,i,alpha,beta)
     3                  +   ym  * EvOpF3(4,i,alpha,beta) )
         enddo
         FKSimulator = norm * FKSimulator
*
****  Dimuon anti-neutrino cross section
*
      elseif(obs(1:9).eq."DIS_DM_NB")then
         call SetProcessDIS("CC")
         call SetProjectileDIS("antineutrino")
         call SetTargetDIS("iron")
         ypc  = yp - 2d0 * ( MN * x * y )**2d0 / Q2
         norm = 100d0 / 2d0 / ( 1d0 + Q2 / MW2 )**2d0 
         do alpha=0,nin(0)
            FKSimulator = FKSimulator + w_int_gen(n,alpha,x)
     1                  * ( ypc * EvOpF2(4,i,alpha,beta) 
     2                  -   y2  * EvOpFL(4,i,alpha,beta)
     3                  -   ym  * EvOpF3(4,i,alpha,beta) )
         enddo
         FKSimulator = norm * FKSimulator
      else
         write(6,*) "In FKSimulator.f:"
         write(6,*) "Invalid observable, obs = ",obs
         write(6,*) "  "
         call exit(-10)
      endif
*
      return
      end
