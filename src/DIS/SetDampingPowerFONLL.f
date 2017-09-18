************************************************************************
*
*     SetDampingPowerFONLL.f:
*
*     This subroutine sets the power with wich the FONLL damping factor
*     suppresses F^{(d)}.
*
************************************************************************
      subroutine SetDampingPowerFONLL(dpc,dpb,dpt)
*
      implicit none
*
      include "../commons/DampingFONLL.h"
*
      integer dpc,dpb,dpt
*
      DampPowerFONLL(4) = dpc
      DampPowerFONLL(5) = dpb
      DampPowerFONLL(6) = dpt
      InDampPowerFONLL  = "done"
*
      return
      end
