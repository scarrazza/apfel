************************************************************************
*
*     SetPolarizedEvolution.f:
*
*     This subroutine enables or disables polarized evolution
*
************************************************************************
      subroutine SetPolarizedEvolution(polev)
*
      implicit none
*
      include "../commons/Polarized.h"
*
      logical polev
*
      Polarized   = polev
      InPolarized = "done"
*
      return
      end
