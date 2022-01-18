************************************************************************
*
*     SetRenFacRatioAlpha.f:
*
*     This subroutine sets the ratio between renormalization and
*     factorization scales in the coupling evolution.
*
************************************************************************
      subroutine SetRenFacRatioAlpha(ratio)
*
      implicit none
*
      include "../commons/krenalpha.h"
*
*     Variables
*
      double precision ratio
*
      krena       = ratio * ratio
      InKrenAlpha = "done"
*
      return
      end
