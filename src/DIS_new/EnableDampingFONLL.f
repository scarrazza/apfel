************************************************************************
*
*     EnableDampingFONLL.f:
*
*     This subroutine enables or disables the damping factor when
*     computing the FONLL structure functions.
*
************************************************************************
      subroutine EnableDampingFONLL(df)
*
      implicit none
*
      include "../commons/DampingFONLL.h"
*
      logical df
*
      DampingFONLL   = df
      InDampingFONLL = "done"
*
      return
      end
