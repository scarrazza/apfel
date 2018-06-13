************************************************************************
*
*     SetAlphaEvolution.f:
*
*     This subroutine sets the solution of the coupling equations.
*
************************************************************************
      subroutine SetAlphaEvolution(ae)
*
      implicit none
*
      include "../commons/AlphaEvolution.h"
*
*     Variables
*
      character*(*) ae
*
      AlphaEvol = trim(ae)
      InAlphaEvol = "done"
*
      return
      end
