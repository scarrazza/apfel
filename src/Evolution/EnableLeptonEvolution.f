************************************************************************
*
*     EnableLeptonEvolution.f:
*
*     This subroutine enables or disables the computation of the global
*     evolution operator that can be used a posteriori to evolve any
*     PDF set.
*
************************************************************************
      subroutine EnableLeptonEvolution(le)
*
      implicit none
*
      include "../commons/LeptEvol.h"
*
      logical le
*
      LeptEvol   = le
      InLeptEvol = "done"
*
      return
      end
