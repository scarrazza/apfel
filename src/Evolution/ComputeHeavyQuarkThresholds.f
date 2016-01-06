************************************************************************
*
*     ComputeHeavyQuarkThresholds.f:
*
*     This subroutine computes the heavy quark thresholds starting from
*     the physical masses and from the ratios bwtween physical masses
*     nad thresholds
*
************************************************************************
      subroutine ComputeHeavyQuarkThresholds
*
      implicit none
*
      include "../commons/m2th.h"
*
*     Variables
*
      integer i
*
      do i=4,6
         m2th(i) = k2th(i) * m2ph(i)
      enddo
*
      return
      end
