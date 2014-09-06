************************************************************************
*
*     SetTheory.f:
*
*     This subroutine sets the FFNS as a default.
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
      character*5 theory
*
      if(theory(1:3).eq."QCD".or.theory(1:3).eq."QED")then
         Th = theory(1:3)
      else
         Th = theory
      endif
      InTheory = "done"
*
      return
      end
