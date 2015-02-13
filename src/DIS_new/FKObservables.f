************************************************************************
*
*     FKObservables.f:
*
*     This function returns simulates the FKgenerator of the NNPDF code.
*
************************************************************************
      function FKObservables(obs,x,Q,y)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/StructureFunctions.h"
      include "../commons/WMass.h"
      include "../commons/ProtonMass.h"
      include "../commons/GFermi.h"
      include "../commons/consts.h"
**
*     Input Variables
*
      double precision x,Q,y
      character*15 obs
**
*     Internal Variables
*
      integer n
      integer alpha
      double precision w_int
      double precision yp,ym,y2,ypc
      double precision Q2,MW2,GF2,MN
      double precision norm
      double precision conv
      parameter(conv=3.893793d10) ! conversion factor from GeV^-2 to 10^-38 cm^2
**
*     Output Variables
*
      double precision FKObservables
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In FKObservables.f:"
         write(6,*) "Invalid value of x =",x
         call exit(-10)
      endif
*
*     Select the grid
*
      do igrid=1,ngrid
         if(x.ge.xmin(igrid).and.x.lt.xmin(igrid+1))then
            goto 101
         endif
      enddo
*
 101  n = inter_degree(igrid)
      FKObservables = 0d0
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
**
*     Electromagnetic Observables
*
****  Light structure function F2light
*
      if(obs(1:7).eq."DIS_F2L")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * F2(3,igrid,alpha)
         enddo
*
****  Up structure function F2u
*
      elseif(obs(1:7).eq."DIS_F2U")then
         write(6,*) "For this observables I don't know yet what to do!"
         call exit(-10)
*
****  Down structure function F2d
*
      elseif(obs(1:7).eq."DIS_F2d")then
         write(6,*) "For this observables I don't know yet what to do!"
         call exit(-10)
*
****  Strange structure function F2s
*
      elseif(obs(1:7).eq."DIS_F2S")then
         write(6,*) "For this observables I don't know yet what to do!"
         call exit(-10)
*
****  Charm structure function F2charm
*
      elseif(obs(1:7).eq."DIS_F2C")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * F2(4,igrid,alpha)
         enddo
*
****  Bottom structure function F2bottom
*
      elseif(obs(1:7).eq."DIS_F2B")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * F2(5,igrid,alpha)
         enddo
*
****  Top structure function F2top
*
      elseif(obs(1:7).eq."DIS_F2T")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * F2(6,igrid,alpha)
         enddo
*
****  Proton structure function F2
*
      elseif(obs(1:7).eq."DIS_F2P")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * F2(7,igrid,alpha)
         enddo
*
****  Deuteron structure function F2
*
      elseif(obs(1:7).eq."DIS_F2D")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * F2(7,igrid,alpha)
         enddo
*
****  Light structure function FLlight
*
      elseif(obs(1:7).eq."DIS_FLL")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * FL(3,igrid,alpha)
         enddo
*
****  Charm structure function FLcharm
*
      elseif(obs(1:7).eq."DIS_FLC")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * FL(4,igrid,alpha)
         enddo
*
****  Bottom structure function FLbottom
*
      elseif(obs(1:7).eq."DIS_FLB")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * FL(5,igrid,alpha)
         enddo
*
****  Top structure function FLtop
*
      elseif(obs(1:7).eq."DIS_FLT")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * FL(6,igrid,alpha)
         enddo
*
****  Proton structure function FL
*
      elseif(obs(1:7).eq."DIS_FLP")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * FL(7,igrid,alpha)
         enddo
*
****  Deuteron structure function FL
*
      elseif(obs(1:7).eq."DIS_FLD")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * FL(7,igrid,alpha)
         enddo
**
*     Neutral-Current Observables
*
****  F2 structure function
*
      elseif(obs(1:10).eq."DIS_F2P_NC")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * F2(7,igrid,alpha)
         enddo
*
****  FL longitudinal structure function
*
      elseif(obs(1:10).eq."DIS_FLP_NC".or.
     1       obs(1:14).eq."DIS_FLP_CON_NC")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * FL(7,igrid,alpha)
         enddo
*
****  F3 structure function
*
      elseif(obs(1:10).eq."DIS_F3P_NC")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * F3(7,igrid,alpha)
         enddo
*
****  Electron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_NCE")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( F2(7,igrid,alpha) 
     2                    - ( y2 / yp ) * FL(7,igrid,alpha)
     3                    + ( ym / yp ) * F3(7,igrid,alpha) )
         enddo
*
****  Positron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_NCP")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( F2(7,igrid,alpha) 
     2                    - ( y2 / yp ) * FL(7,igrid,alpha)
     3                    - ( ym / yp ) * F3(7,igrid,alpha) )
         enddo
*
****  Electron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_NCE_L")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( F2(3,igrid,alpha) 
     2                    - ( y2 / yp ) * FL(3,igrid,alpha)
     3                    + ( ym / yp ) * F3(3,igrid,alpha) )
         enddo
*
****  Positron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_NCP_L")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( F2(3,igrid,alpha) 
     2                    - ( y2 / yp ) * FL(3,igrid,alpha)
     3                    - ( ym / yp ) * F3(3,igrid,alpha) )
         enddo
*
****  Electron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:10).eq."DIS_NCE_CH")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( F2(4,igrid,alpha) 
     2                    - ( y2 / yp ) * FL(4,igrid,alpha)
     3                    + ( ym / yp ) * F3(4,igrid,alpha) )
         enddo
*
****  Positron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:10).eq."DIS_NCP_CH")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( F2(4,igrid,alpha) 
     2                    - ( y2 / yp ) * FL(4,igrid,alpha)
     3                    - ( ym / yp ) * F3(4,igrid,alpha) )
         enddo
*
****  Electron scattering Reduced Cross-Section (bottom)
*
      elseif(obs(1:10).eq."DIS_NCE_BT")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( F2(5,igrid,alpha) 
     2                    - ( y2 / yp ) * FL(5,igrid,alpha)
     3                    + ( ym / yp ) * F3(5,igrid,alpha) )
         enddo
*
****  Positron scattering Reduced Cross-Section (bottom)
*
      elseif(obs(1:10).eq."DIS_NCP_BT")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( F2(5,igrid,alpha) 
     2                    - ( y2 / yp ) * FL(5,igrid,alpha)
     3                    - ( ym / yp ) * F3(5,igrid,alpha) )
         enddo
*
****  Electron scattering Reduced Cross-Section (top)
*
      elseif(obs(1:10).eq."DIS_NCE_TP")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( F2(6,igrid,alpha) 
     2                    - ( y2 / yp ) * FL(6,igrid,alpha)
     3                    + ( ym / yp ) * F3(6,igrid,alpha) )
         enddo
*
****  Positron scattering Reduced Cross-Section (top)
*
      elseif(obs(1:10).eq."DIS_NCP_TP")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( F2(6,igrid,alpha) 
     2                    - ( y2 / yp ) * FL(6,igrid,alpha)
     3                    - ( ym / yp ) * F3(6,igrid,alpha) )
         enddo
*
****  Electron scattering Reduced Cross-Section on deuteron (inclusive)
*
      elseif(obs(1:9).eq."DIS_NCE_D")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( F2(7,igrid,alpha) 
     2                    - ( y2 / yp ) * FL(7,igrid,alpha)
     3                    + ( ym / yp ) * F3(7,igrid,alpha) )
         enddo
*
****  Positron scattering Reduced Cross-Section on deuteron (inclusive)
*
      elseif(obs(1:9).eq."DIS_NCP_D")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( F2(7,igrid,alpha) 
     2                    - ( y2 / yp ) * FL(7,igrid,alpha)
     3                    - ( ym / yp ) * F3(7,igrid,alpha) )
         enddo
**
*     Charged-Current Observables
*
****  Electron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_CCE")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( yp * F2(7,igrid,alpha) 
     2                    -   y2 * FL(7,igrid,alpha)
     3                    -   ym * F3(7,igrid,alpha) )
         enddo
         FKObservables = FKObservables / 4d0
*
****  Positron scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_CCP")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( yp * F2(7,igrid,alpha) 
     2                    -   y2 * FL(7,igrid,alpha)
     3                    +   ym * F3(7,igrid,alpha) )
         enddo
         FKObservables = FKObservables / 4d0
*
****  Electron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_CCE_L")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( yp * F2(3,igrid,alpha) 
     2                    -   y2 * FL(3,igrid,alpha)
     3                    -   ym * F3(3,igrid,alpha) )
         enddo
         FKObservables = FKObservables / 4d0
*
****  Positron scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_CCP_L")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( yp * F2(3,igrid,alpha) 
     2                    -   y2 * FL(3,igrid,alpha)
     3                    +   ym * F3(3,igrid,alpha) )
         enddo
         FKObservables = FKObservables / 4d0
*
****  Electron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_CCE_C")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( yp * F2(4,igrid,alpha) 
     2                    -   y2 * FL(4,igrid,alpha)
     3                    -   ym * F3(4,igrid,alpha) )
         enddo
         FKObservables = FKObservables / 4d0
*
****  Positron scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_CCP_C")then
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( yp * F2(4,igrid,alpha) 
     2                    -   y2 * FL(4,igrid,alpha)
     3                    +   ym * F3(4,igrid,alpha) )
         enddo
         FKObservables = FKObservables / 4d0
**
*     Neutrino Observables
*
****  Neutrino scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_SNU")then
         ypc  = yp - 2d0 * ( MN * x * y )**2d0 / Q2
         norm = GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2d0 )
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( ypc * F2(7,igrid,alpha) 
     2                    -   y2  * FL(7,igrid,alpha)
     3                    +   ym  * F3(7,igrid,alpha) )
         enddo
         FKObservables = norm * FKObservables
*
****  Antineutrino scattering Reduced Cross-Section (inclusive)
*
      elseif(obs(1:7).eq."DIS_SNB")then
         ypc  = yp - 2d0 * ( MN * x * y )**2d0 / Q2
         norm = GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2d0 )
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( ypc * F2(7,igrid,alpha) 
     2                    -   y2  * FL(7,igrid,alpha)
     3                    -   ym  * F3(7,igrid,alpha) )
         enddo
         FKObservables = norm * FKObservables
*
****  Neutrino scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_SNU_L")then
         ypc  = yp - 2d0 * ( MN * x * y )**2d0 / Q2
         norm = GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2d0 )
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( ypc * F2(3,igrid,alpha) 
     2                    -   y2  * FL(3,igrid,alpha)
     3                    +   ym  * F3(3,igrid,alpha) )
         enddo
         FKObservables = norm * FKObservables
*
****  Antineutrino scattering Reduced Cross-Section (light)
*
      elseif(obs(1:9).eq."DIS_SNB_L")then
         ypc  = yp - 2d0 * ( MN * x * y )**2d0 / Q2
         norm = GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2d0 )
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( ypc * F2(3,igrid,alpha) 
     2                    -   y2  * FL(3,igrid,alpha)
     3                    -   ym  * F3(3,igrid,alpha) )
         enddo
         FKObservables = norm * FKObservables
*
****  Neutrino scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_SNU_C")then
         ypc  = yp - 2d0 * ( MN * x * y )**2d0 / Q2
         norm = GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2d0 )
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( ypc * F2(4,igrid,alpha) 
     2                    -   y2  * FL(4,igrid,alpha)
     3                    +   ym  * F3(4,igrid,alpha) )
         enddo
         FKObservables = norm * FKObservables
*
****  Antineutrino scattering Reduced Cross-Section (charm)
*
      elseif(obs(1:9).eq."DIS_SNB_C")then
         ypc  = yp - 2d0 * ( MN * x * y )**2d0 / Q2
         norm = GF2 * MN / ( 2d0 * pi * ( 1d0 + Q2 / MW2 )**2d0 )
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( ypc * F2(4,igrid,alpha) 
     2                    -   y2  * FL(4,igrid,alpha)
     3                    -   ym  * F3(4,igrid,alpha) )
         enddo
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
         ypc  = yp - 2d0 * ( MN * x * y )**2d0 / Q2
         norm = 100d0 / 2d0 / ( 1d0 + Q2 / MW2 )**2d0 
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( ypc * F2(4,igrid,alpha) 
     2                    -   y2  * FL(4,igrid,alpha)
     3                    +   ym  * F3(4,igrid,alpha) )
         enddo
         FKObservables = norm * FKObservables
*
****  Dimuon anti-neutrino cross section
*
      elseif(obs(1:9).eq."DIS_DM_NB")then
         ypc  = yp - 2d0 * ( MN * x * y )**2d0 / Q2
         norm = 100d0 / 2d0 / ( 1d0 + Q2 / MW2 )**2d0 
         do alpha=0,nin(0)
            FKObservables = FKObservables + w_int(n,alpha,x)
     1                    * ( ypc * F2(4,igrid,alpha) 
     2                    -   y2  * FL(4,igrid,alpha)
     3                    -   ym  * F3(4,igrid,alpha) )
         enddo
         FKObservables = norm * FKObservables
      else
         write(6,*) "In FKObservables.f:"
         write(6,*) "Invalid observable, obs = ",obs
         write(6,*) "  "
         call exit(-10)
      endif
*
      return
      end
