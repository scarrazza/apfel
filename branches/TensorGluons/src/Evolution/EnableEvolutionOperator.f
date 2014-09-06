************************************************************************
*
*     EnableEvolutionOperator.f:
*
*     This subroutine enables or disables the computation of the global
*     evolution operator that can be used a posteriori to evolve any
*     PDF set.
*
************************************************************************
      subroutine EnableEvolutionOperator(eo)
*
      implicit none
*
      include "../commons/EvolOp.h"
*
      logical eo
*
      EvolOp   = eo
      InEvolOp = "done"
*
      return
      end
