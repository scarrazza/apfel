************************************************************************
*
*     ComputeChargesDIS.f:
*
*     It sets the charges used for the computation of the DIS observables.
*
************************************************************************
      subroutine ComputeChargesDIS(Q2,bq,dq,bqt)
*
      implicit none
*
      include "../commons/ProcessDIS.h"
      include "../commons/PolarizationDIS.h"
      include "../commons/ProjectileDIS.h"
      include "../commons/Sin2ThetaW.h"
      include "../commons/ZedMass.h"
      include "../commons/SelectedCharge.h"
      include "../commons/TimeLike.h"
      include "../commons/m2th.h"
      include "../commons/MaxFlavourPDFs.h"
      include "../commons/PropagatorCorrection.h"
      include "../commons/EWCouplings.h"
      include "../commons/NCComponent.h"
**
*     Input Variables
*
      double precision Q2
**
*     Input Variables
*
      integer i
      integer ie
      integer nfi,nff,nf
      double precision pol
      double precision eq(6),eq2(6)
      double precision vq(6),aq(6)
      double precision ve,ae
      double precision pz,pz2
      double precision GammaZ
      double precision sumbq,sumdq,sumbqt
      parameter(GammaZ = 2.4952d0)
**
*     Double precision
*
      double precision bq(0:6),dq(0:6),bqt(0:6)
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
      eq(3) = eq(1) ! - 1d0 / 3d0
      eq(4) = eq(2) ! 2d0 / 3d0
      eq(5) = eq(1) ! - 1d0 / 3d0
      eq(6) = eq(2) ! 2d0 / 3d0
*
*     Squared Charges
*
      eq2(1) = eq(1) * eq(1)
      eq2(2) = eq(2) * eq(2)
      eq2(3) = eq2(1) ! eq(3) * eq(3)
      eq2(4) = eq2(2) ! eq(4) * eq(4)
      eq2(5) = eq2(1) ! eq(5) * eq(5)
      eq2(6) = eq2(2) ! eq(6) * eq(6)
*
*     Vector Couplings
*
      if(ExtCoup)then
         vq(1) = VectorD
         vq(2) = VectorU
      else
         vq(1) = - 0.5d0 + 2d0 / 3d0 * Sin2ThetaW
         vq(2) = + 0.5d0 - 4d0 / 3d0 * Sin2ThetaW
      endif
      vq(3) = vq(1) ! - 0.5d0 + 2d0 / 3d0 * Sin2ThetaW
      vq(4) = vq(2) ! + 0.5d0 - 4d0 / 3d0 * Sin2ThetaW
      vq(5) = vq(1) ! - 0.5d0 + 2d0 / 3d0 * Sin2ThetaW
      vq(6) = vq(2) ! + 0.5d0 - 4d0 / 3d0 * Sin2ThetaW
*
*     Axial Couplings
*
      if(ExtCoup)then
         aq(1) = AxialD
         aq(2) = AxialU
      else
         aq(1) = - 0.5d0
         aq(2) = + 0.5d0
      endif
      aq(3) = aq(1) ! - 0.5d0
      aq(4) = aq(2) ! + 0.5d0
      aq(5) = aq(1) ! - 0.5d0
      aq(6) = aq(2) ! + 0.5d0
*
*     Vector and Axial Electron Couplings
*
      ve = - 0.5d0 + 2d0 * Sin2ThetaW
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
         ve = 0.5d0
         ae = 0.5d0
      endif
*
*     Initialize charges
*
      do i=0,6
         bq(i)  = 0d0
         dq(i)  = 0d0
         bqt(i) = 0d0
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
            bq(i)  = eq2(i)
            dq(i)  = 0d0
            bqt(i) = eq2(i)
         enddo
      elseif(ProcessDIS.eq."NC")then
         if(TimeLike)then
            pz  = Q2 * ( Q2 -  MZ**2 )
     1          / ( ( Q2 - MZ**2 )**2 + ( MZ * GammaZ )**2 ) 
     2          / ( 4d0 * Sin2ThetaW * ( 1d0 - Sin2ThetaW ) )
            pz2 = Q2**2 
     1          / ( ( Q2 - MZ**2 )**2 + ( MZ * GammaZ )**2 ) 
     2          / ( ( 4d0 * Sin2ThetaW * ( 1d0 - Sin2ThetaW ) ) )**2
         else
            pz  = Q2 / ( Q2 + MZ**2 ) 
     1          / ( 4d0 * Sin2ThetaW * ( 1d0 - Sin2ThetaW ) )
            pz2 = pz * pz
         endif
*     Apply propagator correction
         pz  = pz  / ( 1d0 - DeltaR )
         pz2 = pz2 / ( 1d0 - DeltaR )**2
         if(NCComponent.eq."gg")then
            do i=nfi,nff
               bq(i)  = eq2(i)
               dq(i)  = 0d0
               bqt(i) = eq2(i)
            enddo
         elseif(NCComponent.eq."gZ")then
            do i=nfi,nff
c               bq(i)  = - 2d0*eq(i)*vq(i) * ( ve + ie * pol * ae ) * pz
c               dq(i)  = - 2d0*eq(i)*aq(i) * ( ae + ie * pol * ve ) * pz
c               bqt(i) = - 2d0*eq(i)*vq(i) * ( ve + ie * pol * ae ) * pz
               bq(i)  = 2d0 * eq(i) * vq(i)
               dq(i)  = 2d0 * eq(i) * aq(i)
               bqt(i) = 2d0 * eq(i) * vq(i)
            enddo
         elseif(NCComponent.eq."ZZ")then
            do i=nfi,nff
c               bq(i)  = ( ve**2 + ae**2 + ie * pol * 2d0 * ve * ae )
c     1              * ( vq(i)**2 + aq(i)**2 ) * pz2
c               dq(i)  = 2d0 * vq(i) * aq(i) * ( 2d0 * ve * ae
c     1              + ie * pol * ( ve**2 + ae**2 ) ) * pz2
c               bqt(i) = ( ve**2 + ae**2 + ie * pol * 2d0 * ve * ae )
c     1              * ( vq(i)**2 - aq(i)**2 ) * pz2
               bq(i)  = vq(i)**2 + aq(i)**2
               dq(i)  = 2d0 * vq(i) * aq(i)
               bqt(i) = vq(i)**2 - aq(i)**2
            enddo
         else
            do i=nfi,nff
               bq(i)  = eq2(i)
     1              - 2d0 * eq(i) * vq(i) * ( ve + ie * pol * ae ) * pz
     2              + ( ve**2 + ae**2 + ie * pol * 2d0 * ve * ae )
     3              * ( vq(i)**2 + aq(i)**2 ) * pz2
               dq(i)  = - 2d0*eq(i)*aq(i) * ( ae + ie * pol * ve ) * pz
     1              + 2d0 * vq(i) * aq(i) * ( 2d0 * ve * ae
     2              + ie * pol * ( ve**2 + ae**2 ) ) * pz2
               bqt(i) = eq2(i)
     1              - 2d0 * eq(i) * vq(i) * ( ve + ie * pol * ae ) * pz
     2              + ( ve**2 + ae**2 + ie * pol * 2d0 * ve * ae )
     3              * ( vq(i)**2 - aq(i)**2 ) * pz2
            enddo
         endif
      endif
*
*     Normalize charges for the SIA structure functions
*
      if(TimeLike)then
         if(Q2.ge.m2th(6))then
            nf = 6
         elseif(Q2.ge.m2th(5))then
            nf = 5
         elseif(Q2.ge.m2th(4))then
            nf = 4
         else
            nf = 3
         endif
         if(nf.gt.nfMaxPDFs) nf = nfMaxPDFs
*
         sumbq  = 0d0
         sumdq  = 0d0
         sumbqt = 0d0
         do i=1,nf
            sumbq  = sumbq + bq(i)
            sumdq  = sumdq + dq(i)
            sumbqt = sumdq + bqt(i)
         enddo
*
         bq(0)  = sumbq
         dq(0)  = sumdq
         bqt(0) = sumbqt
         do i=1,6
            bq(i)  = bq(i)  / sumbq
            dq(i)  = dq(i)  / sumdq
            bqt(i) = bqt(i) / sumbqt
         enddo
      endif
*
      return
      end
