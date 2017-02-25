************************************************************************
*
*     EnableSFNLOQEDCorrections.f:
*
*     This subroutine enables or disables the O(alpha) corrections
*     in the DIS sturcture functions. This is active only if the "QUniD"
*     theory is used.
*
************************************************************************
      subroutine EnableSFNLOQEDCorrections(qedsfc)
*
      implicit none
*
      include "../commons/NLOQEDCorrections.h"
*
      logical qedsfc
*
      SFNLOQED   = qedsfc
      InSFNLOQED = "done"
*
      return
      end
