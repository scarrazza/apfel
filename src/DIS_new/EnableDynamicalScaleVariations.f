************************************************************************
*
*     EnableDynamicalScaleVariations.f
*
*     This subroutine enables or disables the possibility to perform
*     factorization and renormalization scale variations point by point
*     without requiring the ratio \mu_{R,F} / Q to be constant.
*     Limitations: \mu_F = \mu_R and the code slows down as the additional
*     cotrubutions to the coefficient functions are included in real time.
*
************************************************************************
      subroutine EnableDynamicalScaleVariations(dsv)
*
      implicit none
*
      include "../commons/DynScVar.h"
*
      logical dsv
*
      DynScVar   = dsv
      InDynScVar = "done"
*
      return
      end
