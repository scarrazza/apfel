************************************************************************
*
*     SetSinThetaW.f:
*
*     This subroutine sets sin(\theta_W).
*
************************************************************************
      subroutine SetSinThetaW(sw)
*
      implicit none
*
      include "../commons/SinThetaW.h"
*
*     Variables
*
      double precision sw
*
      SinThetaW   = sw
      InSinThetaW = "done"
*
      return
      end
