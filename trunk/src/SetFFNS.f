************************************************************************
*
*     SetFFNS.f:
*
*     This subroutine sets the FFNS as a default.
*
************************************************************************
      subroutine SetFFNS(nfl)
*
      implicit none
*
      include "../commons/Evs.h"
      include "../commons/Nf_FF.h"
*
*     Variables
*
      integer nfl
*
      Evs   = "FF"
      Nf_FF = nfl
      InEvs = "done"
*
      return
      end
