************************************************************************
*
*     initDYcouplings.f:
*
*     Initialize the Drell-Yan couplings.
*
************************************************************************
      subroutine initDYcouplings
*
      implicit none
*
      include "../commons/DYcouplings.h"
**
*     Internal Variables
*
      integer ibos,iq,jq
      integer nff
      double precision sw,GetSin2ThetaW
      double precision V_UD,V_US,V_CD,V_CS,GetCKM
      double precision V_CKM(-6:6,-6:6)
*
*     sin^2(\theta_W)
*
      sw = GetSin2ThetaW()
*
      do ibos=1,4
         do iq=-6,6
            VV(iq,ibos) = 0d0
            AA(iq,ibos) = 0d0
            do jq=-6,6
               CII(iq,jq,ibos) = 0d0
               CIF(iq,jq,ibos) = 0d0
               CFF(iq,jq,ibos) = 0d0
            enddo
         enddo
      enddo   
*
      do iq=-6,6
         do jq=-6,6
            V_CKM(iq,jq) = 0d0
         enddo
      enddo
*
*     CKM matrix elements
*
      V_UD = GetCKM(1,1)
      V_US = GetCKM(1,2)
      V_CD = GetCKM(2,1)
      V_CS = GetCKM(2,2)
*     Vud
      V_CKM(-1,2) = V_UD
      V_CKM(-2,1) = V_UD
      V_CKM(1,-2) = V_UD
      V_CKM(2,-1) = V_UD
      V_CKM(1,2)  = V_UD
      V_CKM(2,1)  = V_UD
      V_CKM(-1,-2)= V_UD
      V_CKM(-2,-1)= V_UD
*     Vus
      V_CKM(-2,3) = V_US
      V_CKM(-3,2) = V_US
      V_CKM(2,-3) = V_US
      V_CKM(3,-2) = V_US
      V_CKM(2,3)  = V_US
      V_CKM(3,2)  = V_US
      V_CKM(-2,-3)= V_US
      V_CKM(-3,-2)= V_US
*     Vcd
      V_CKM(-1,4) = V_CD
      V_CKM(-4,1) = V_CD
      V_CKM(1,-4) = V_CD
      V_CKM(4,-1) = V_CD
      V_CKM(1,4)  = V_CD
      V_CKM(4,1)  = V_CD
      V_CKM(-1,-4)= V_CD
      V_CKM(-4,-1)= V_CD
*     Vcs
      V_CKM(-3,4) = V_CS
      V_CKM(-4,3) = V_CS
      V_CKM(3,-4) = V_CS
      V_CKM(4,-3) = V_CS
      V_CKM(3,4)  = V_CS
      V_CKM(4,3)  = V_CS
      V_CKM(-3,-4)= V_CS
      V_CKM(-4,-3)= V_CS
*     
*     Electroweak couplings
*
*     ibos = 1 --> Gamma
*
      VV(1,1) = - 1d0 / 3d0        !down
      VV(2,1) =   2d0 / 3d0        !up
      VV(3,1) = - 1d0 / 3d0        !strange
      VV(4,1) =   2d0 / 3d0        !charm
      VV(5,1) = - 1d0 / 3d0        !bottom
      VV(6,1) =   2d0 / 3d0        !top
*
*     ibos = 2 --> Z
*
      VV(1,2) = ( - 1d0 + 4d0 / 3d0 * SW ) / 2d0
      VV(2,2) = (   1d0 - 8d0 / 3d0 * SW ) / 2d0
      VV(3,2) = ( - 1d0 + 4d0 / 3d0 * SW ) / 2d0
      VV(4,2) = (   1d0 - 8d0 / 3d0 * SW ) / 2d0
      VV(5,2) = ( - 1d0 + 4d0 / 3d0 * SW ) / 2d0
      VV(6,2) = (   1d0 - 8d0 / 3d0 * SW ) / 2d0
*
      AA(1,2)  =   1d0 / 2d0
      AA(2,2)  = - 1d0 / 2d0
      AA(3,2)  =   1d0 / 2d0
      AA(4,2)  = - 1d0 / 2d0
      AA(5,2)  =   1d0 / 2d0
      AA(6,2)  = - 1d0 / 2d0
*
*     ibos = 3 --> W+
*
      VV(1,3)  = 1d0 / DSQRT(2d0)
      VV(2,3)  = 1d0 / DSQRT(2d0)
      VV(3,3)  = 1d0 / DSQRT(2d0)
      VV(4,3)  = 1d0 / DSQRT(2d0)
      VV(5,3)  = 1d0 / DSQRT(2d0)
      VV(6,3)  = 1d0 / DSQRT(2d0)
*
      AA(1,3)  = - 1d0 / DSQRT(2d0)
      AA(2,3)  = - 1d0 / DSQRT(2d0)
      AA(3,3)  = - 1d0 / DSQRT(2d0)
      AA(4,3)  = - 1d0 / DSQRT(2d0)
      AA(5,3)  = - 1d0 / DSQRT(2d0)
      AA(6,3)  = - 1d0 / DSQRT(2d0)
*
*     ibos = 4 --> W-
*
      VV(1,4)  = 1d0 / DSQRT(2d0)
      VV(2,4)  = 1d0 / DSQRT(2d0)
      VV(3,4)  = 1d0 / DSQRT(2d0)
      VV(4,4)  = 1d0 / DSQRT(2d0)
      VV(5,4)  = 1d0 / DSQRT(2d0)
      VV(6,4)  = 1d0 / DSQRT(2d0)
*
      AA(1,4)  = - 1d0 / DSQRT(2d0)
      AA(2,4)  = - 1d0 / DSQRT(2d0)
      AA(3,4)  = - 1d0 / DSQRT(2d0)
      AA(4,4)  = - 1d0 / DSQRT(2d0)
      AA(5,4)  = - 1d0 / DSQRT(2d0)
      AA(6,4)  = - 1d0 / DSQRT(2d0)
*
*     All the couplings get a minus under charge conjugation
*
      do ibos=1,4
         do jq=1,6
            VV(-jq,ibos) = - VV(jq,ibos)
            AA(-jq,ibos) = - AA(jq,ibos)
         enddo
      enddo
*
*     Fill the Drell-Yan C matrices (in the notation of Anastasiou et al.)
*
      do iq=-6,6
         CII(iq,-iq,1) = 1d0
         CII(iq,-iq,2) = 1d0
         CFF(iq,-iq,1) = 1d0
         CFF(iq,-iq,2) = 1d0
         CIF(iq,iq,1)  = 1d0
         CIF(iq,iq,2)  = 1d0
         do jq=-6,6
            if(DABS(VV(iq,1)+VV(jq,1)-1d0).lt.1d-8)then 
               CII(iq,jq,3) = V_CKM(iq,jq)**2d0
               CFF(iq,jq,4) = V_CKM(iq,jq)**2d0
            endif
            if(DABS(VV(iq,1)+VV(jq,1)+1d0).lt.1d-8)then 
               CII(iq,jq,4) = V_CKM(iq,jq)**2d0
               CFF(iq,jq,3) = V_CKM(iq,jq)**2d0
            endif
            if(DABS(VV(iq,1)-VV(jq,1)-1d0).lt.1d-8)then 
               CIF(iq,jq,3) = V_CKM(iq,jq)**2d0
            endif
            if(DABS(VV(iq,1)-VV(jq,1)+1d0).lt.1d-8)then 
               CIF(iq,jq,4) = V_CKM(iq,jq)**2d0
            endif
         enddo
      enddo
*
      do ibos=1,4
         do iq=-6,6
            CII(0,iq,ibos) = 0d0
            CII(iq,0,ibos) = 0d0
            CIF(0,iq,ibos) = 0d0
            CIF(iq,0,ibos) = 0d0
            CFF(0,iq,ibos) = 0d0
            CFF(iq,0,ibos) = 0d0
         enddo
      enddo
*
*     Compute simplified combinations for NLO only
*
      do ibos=1,4
         do nff=1,6
            do iq=-nff,nff
               CIF_NLO(nff,iq,ibos) = 0d0
               do jq = -nff,nff
                  CIF_NLO(nff,iq,ibos) = CIF_NLO(nff,iq,ibos) 
     1                                 + CIF(iq,jq,ibos)
               enddo
            enddo
         enddo
      enddo
*
      return 
      end
