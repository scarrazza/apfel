************************************************************************
*
*     initHELL.f:
*
*     This routine initializes the HELL code by Bonvini for the small-x
*     resummation.
*
************************************************************************
      subroutine initHELL
*
      implicit none
*
      include "../commons/Smallx.h"
      include "../commons/ipt.h"
*
      call HELLLogOrder(LogAcc)
      call HELLOrder(ipt)
      call HELL
*
      return
      end
