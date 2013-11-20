************************************************************************
*
*     SetRenFacRatio.f:
*
*     This subroutine sets the ratio between renormalization and
*     factiorization scales.
*
************************************************************************
      subroutine SetRenFacRatio(ratio)
*
      implicit none
*
      include "../commons/kren.h"
*
*     Variables
*
      double precision ratio
*
      kren   = ratio * ratio
      InKren = "done"
*
      return
      end
