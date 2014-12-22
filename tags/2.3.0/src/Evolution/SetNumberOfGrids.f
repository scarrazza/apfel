************************************************************************
*
*     SetNumberOfGrids.f:
*
*     This subroutine sets the number of x-space grids that will be used
*     in the computation
*
************************************************************************
      subroutine SetNumberOfGrids(n)
*
      implicit none
*
      include "../commons/grid.h"
*
*     Variables
*
      integer n
*
      ThereAreExtGrids = .false.
*
      ngrid  = n
      InGrid = "done"
*
      return
      end
