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
      include "../commons/SinThetaW.h"
      include "../commons/ZedMass.h"
**
*     Input Variables
*
      double precision Q2
**
*     Input Variables
*
      integer i
      integer ie
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
*     temporary definitions
*
      ie  = - 1   ! Electron
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
      eq2(1) = 1d0 / 9d0
      eq2(2) = 4d0 / 9d0
      eq2(3) = 1d0 / 9d0
      eq2(4) = 4d0 / 9d0
      eq2(5) = 1d0 / 9d0
      eq2(6) = 4d0 / 9d0
*
*     Vector Couplings
*
      vq(1) = - 0.5d0 + 2d0 / 3d0 * sw
      vq(2) = + 0.5d0 - 4d0 / 3d0 * sw
      vq(3) = - 0.5d0 + 2d0 / 3d0 * sw
      vq(4) = + 0.5d0 - 4d0 / 3d0 * sw
      vq(5) = - 0.5d0 + 2d0 / 3d0 * sw
      vq(6) = + 0.5d0 - 4d0 / 3d0 * sw
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
      ve = - 0.5d0 + 2d0 * sw
      ae = - 0.5d0
*
      if(ProcessDIS.eq."EM")then
         do i=1,6
            bq(i) = eq2(i)
            dq(i) = 0d0
         enddo
      elseif(ProcessDIS.eq."NC")then
         pz = Q2 / ( Q2 + MZ**2d0 ) / ( 4d0 * sw * ( 1d0 - sw ) )
         do i=1,6
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
