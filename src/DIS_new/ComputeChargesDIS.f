************************************************************************
*
*     ComputeChargesDIS.f:
*
*     It sets the charges used for the computation of the DIS observables.
*
************************************************************************
      subroutine ComputeChargesDIS(Q2,bq,dq)
*
      implicit none
*
      include "../commons/ProcessDIS.h"
      include "../commons/PolarizationDIS.h"
      include "../commons/ProjectileDIS.h"
      include "../commons/SinThetaW.h"
      include "../commons/ZedMass.h"
      include "../commons/SelectedCharge.h"
**
*     Input Variables
*
      double precision Q2
**
*     Input Variables
*
      integer i
      integer ie
      integer nfi,nff
      double precision pol
      double precision eq(6),eq2(6)
      double precision vq(6),aq(6)
      double precision ve,ae
      double precision pz
**
*     Double precision
*
      double precision bq(6),dq(6)
*
*     Polarization fraction
*
      pol = PolarizationDIS
*
*     Projectile
*
      if(ProjectileDIS(1:8).eq."electron")then
         ie  = - 1
      elseif(ProjectileDIS(1:8).eq."positron")then
         ie  = 1
      elseif(ProjectileDIS(1:8).eq."neutrino")then
         ie  = 1
      elseif(ProjectileDIS.eq."antineutrino")then
         ie  = - 1
      endif
*
*     Initialize charges and couplings
*
*     Electric Charges
*
      eq(1) = - 1d0 / 3d0
      eq(2) = 2d0 / 3d0
      eq(3) = - 1d0 / 3d0
      eq(4) = 2d0 / 3d0
      eq(5) = - 1d0 / 3d0
      eq(6) = 2d0 / 3d0
*
*     Squared Charges
*
      eq2(1) = eq(1) * eq(1) ! 1d0 / 9d0
      eq2(2) = eq(2) * eq(2) ! 4d0 / 9d0
      eq2(3) = eq(3) * eq(3) ! 1d0 / 9d0
      eq2(4) = eq(4) * eq(4) ! 4d0 / 9d0
      eq2(5) = eq(5) * eq(5) ! 1d0 / 9d0
      eq2(6) = eq(6) * eq(6) ! 4d0 / 9d0
*
*     Vector Couplings
*
      vq(1) = - 0.5d0 + 2d0 / 3d0 * SinThetaW
      vq(2) = + 0.5d0 - 4d0 / 3d0 * SinThetaW
      vq(3) = - 0.5d0 + 2d0 / 3d0 * SinThetaW
      vq(4) = + 0.5d0 - 4d0 / 3d0 * SinThetaW
      vq(5) = - 0.5d0 + 2d0 / 3d0 * SinThetaW
      vq(6) = + 0.5d0 - 4d0 / 3d0 * SinThetaW
*
*     Axial Couplings
*
      aq(1) = - 0.5d0
      aq(2) = + 0.5d0
      aq(3) = - 0.5d0
      aq(4) = + 0.5d0
      aq(5) = - 0.5d0
      aq(6) = + 0.5d0
*
*     Vector and Axial Electron Couplings
*
      ve = - 0.5d0 + 2d0 * SinThetaW
      ae = - 0.5d0
*
*     Set the electric charges to zero if the projectile is either
*     a neutrino or an antineutrino and correct the vector and axial
*     couplings.
*
      if(ProjectileDIS(1:8).eq."neutrino".or.
     1   ProjectileDIS.eq."antineutrino")then
         do i=1,6
            eq(i)  = 0d0
            eq2(i) = 0d0
         enddo
         ve = 0.5d0 + 2d0 * SinThetaW
         ae = 0.5d0
      endif
*
*     Initialize charges
*
      do i=1,6
         bq(i) = 0d0
         dq(i) = 0d0
      enddo
*
*     Select the charge
*
      nfi = 1
      nff = 6
      if(SelectedCharge(1:4).eq."down")then
         nfi = 1
         nff = 1
      elseif(SelectedCharge(1:2).eq."up")then
         nfi = 2
         nff = 2
      elseif(SelectedCharge(1:7).eq."strange")then
         nfi = 3
         nff = 3
      elseif(SelectedCharge(1:5).eq."charm")then
         nfi = 4
         nff = 4
      elseif(SelectedCharge(1:6).eq."bottom")then
         nfi = 5
         nff = 5
      elseif(SelectedCharge(1:3).eq."top")then
         nfi = 6
         nff = 6
      endif
*
      if(ProcessDIS.eq."EM")then
         do i=nfi,nff
            bq(i) = eq2(i)
            dq(i) = 0d0
         enddo
      elseif(ProcessDIS.eq."NC")then
         pz = Q2 / ( Q2 + MZ**2d0 ) 
     1      / ( 4d0 * SinThetaW * ( 1d0 - SinThetaW ) )
         do i=nfi,nff
            bq(i) = eq2(i) 
     1            - 2d0 * eq(i) * vq(i) * ( ve + ie * pol * ae ) * pz
     2            + ( ve**2d0 + ae**2d0 ) * ( vq(i)**2d0 + aq(i)**2d0 
     3            + ie * pol * 2d0 * ve * ae ) * pz**2d0 
            dq(i) = - 2d0 * eq(i) * aq(i) * ( ae + ie * pol * ve ) * pz
     1            + 2d0 * vq(i) * aq(i) * ( 2d0 * ve * ae 
     2            + ie * pol * ( ve**2d0 + ae**2d0 ) ) * pz**2d0
         enddo
      endif
*
      return
      end
