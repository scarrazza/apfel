************************************************************************
*
*     EnableTargetMassCorrections.f:
*
*     This subroutine enables or disables the computation of the target
*     mass corrections due to the finite mass of the proton for the DIS
*     structure functions.
*
************************************************************************
      subroutine EnableTargetMassCorrections(tc)
*
      implicit none
*
      include "../commons/TMC.h"
*
      logical tc
*
      TMC   = tc
      InTMC = "done"
*
      return
      end
