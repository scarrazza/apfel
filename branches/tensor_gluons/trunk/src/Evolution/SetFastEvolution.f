************************************************************************
*
*     SetFastEvolution.f:
*
*     This subroutine enables or disables the fast evolution.
*
************************************************************************
      subroutine SetFastEvolution(fe)
*
      implicit none
*
      include "../commons/FastEvol.h"
*
      logical fe
*
      FastEvol   = fe
      InFastEvol = "done"
*
      return
      end
