************************************************************************
*
*     SelectCharge.f:
*
*     This subroutine selects one particular charge in the NC structure
*     functions (needed fo selecting up down and strange components).
*
************************************************************************
      subroutine SelectCharge(selch)
*
      implicit none
*
      include "../commons/SelectedCharge.h"
*
*     Variables
*
      character*(*) selch
*
      SelectedCharge = trim(selch)
      InSelectedCharge = "done"
*
      return
      end
