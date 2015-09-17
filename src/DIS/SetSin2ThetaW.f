************************************************************************
*
*     SetSin2ThetaW.f:
*
*     This subroutine sets sin^2(\theta_W).
*
************************************************************************
      subroutine SetSin2ThetaW(sw)
*
      implicit none
*
      include "../commons/Sin2ThetaW.h"
*
*     Variables
*
      double precision sw
*
      Sin2ThetaW   = sw
      InSin2ThetaW = "done"
*
      return
      end
