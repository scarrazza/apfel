************************************************************************
*
*     GetSinThetaW.f:
*
*     This function returns sin(\theta_W).
*
************************************************************************
      function GetSinThetaW()
*
      implicit none
*
      include "../commons/SinThetaW.h"
*
*     Variables
*
      double precision GetSinThetaW
*
      GetSinThetaW = SinThetaW
*
      return
      end
