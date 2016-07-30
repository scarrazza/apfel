************************************************************************
*
*     EnableNLOQEDCorrections.f:
*
*     This subroutine enables or disables the O(alpha_s alpha) and
*     O(alpha^2) corrections. This is active only if the "QUniD" theory
*     is used.
*
************************************************************************
      subroutine EnableNLOQEDCorrections(qedc)
*
      implicit none
*
      include "../commons/NLOQEDCorrections.h"
*
      logical qedc
*
      NLOQED   = qedc
      InNLOQED = "done"
*
      return
      end
