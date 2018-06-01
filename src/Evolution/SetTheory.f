************************************************************************
*
*     SetTheory.f:
*
*     This subroutine sets the evolution theory to be used.
*
************************************************************************
      subroutine SetTheory(theory)
*
      implicit none
*
      include "../commons/Th.h"
*
*     Variables
*
      character*(*) theory
      Th = trim(theory)
      InTheory = "done"
*
      return
      end
