************************************************************************
*
*     SetVFNS.f:
*
*     This subroutine sets the VFNS as a default.
*
************************************************************************
      subroutine SetVFNS
*
      implicit none
*
      include "../commons/Evs.h"
*
      Evs   = "VF"
      InEvs = "done"
*
      return
      end
