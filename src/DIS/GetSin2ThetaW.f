************************************************************************
*
*     GetSin2ThetaW.f:
*
*     This function returns sin^2(\theta_W).
*
************************************************************************
      function GetSin2ThetaW()
*
      implicit none
*
      include "../commons/Sin2ThetaW.h"
*
*     Variables
*
      double precision GetSin2ThetaW
*
      if(InSin2ThetaW.ne."done")then
         write(6,*) "GetSin2ThetaW: Parameter not initialized"
         write(6,*) "Set it by means of 'SetSin2ThetaW'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      GetSin2ThetaW = Sin2ThetaW
*
      return
      end
