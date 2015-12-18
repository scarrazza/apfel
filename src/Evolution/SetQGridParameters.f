************************************************************************
*
*     SetQGridParameters.f:
*
*     This subroutine sets the the parameters of the Q2-space grid.
*
************************************************************************
      subroutine SetQGridParameters(npQ,degQ)
*
      implicit none
*
      include "../commons/gridQ.h"
*
*     Variables
*
      integer npQ,degQ
*
      nQ2g          = npQ
      inter_degreeQ = degQ
      InQGrid       = "done"
*
      return
      end
