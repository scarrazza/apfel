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
      if(inSinThetaW.ne."done")then
         write(6,*) "GetSinThetaW: Parameter not initialized"
         write(6,*) "Set it by means of 'SetSinThetaW'"
         write(6,*) "  "
         call exit(-10)
      endif
*
      GetSinThetaW = SinThetaW
*
      return
      end
