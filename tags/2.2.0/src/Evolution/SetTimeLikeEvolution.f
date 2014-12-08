************************************************************************
*
*     SetTimeLikeEvolution.f:
*
*     This subroutine enables or disables time-like evolution for the
*     fragmentation functions.
*
************************************************************************
      subroutine SetTimeLikeEvolution(tl)
*
      implicit none
*
      include "../commons/TimeLike.h"
*
      logical tl
*
      TimeLike   = tl
      InTimeLike = "done"
*
      return
      end
