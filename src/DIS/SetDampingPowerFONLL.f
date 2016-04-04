************************************************************************
*
*     SetDampingPowerFONLL.f:
*
*     This subroutine sets the power with wich the FONLL damping factor
*     suppresses F^{(d)}.
*
************************************************************************
      subroutine SetDampingPowerFONLL(dp)
*
      implicit none
*
      include "../commons/DampingFONLL.h"
*
      integer dp
*
      DampPowerFONLL   = dp
      InDampPowerFONLL = "done"
*
      return
      end
