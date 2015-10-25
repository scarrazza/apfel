************************************************************************
*
*     EnableIntrinsicCharm.f:
*
*     This subroutine enables or disables the computation of the
*     intrinsic charm component to the massive DIS structure functions.
*
************************************************************************
      subroutine EnableIntrinsicCharm(ic)
*
      implicit none
*
      include "../commons/IntrinsicCharm.h"
*
      logical ic
*
      IntrinsicCharm   = ic
      InIntrinsicCharm = "done"
*
      return
      end
