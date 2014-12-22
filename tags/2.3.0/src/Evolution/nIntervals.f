************************************************************************
*
*     nIntervals.f:
*
*     This function returns the number nodes of the joint x-space grid.
*
************************************************************************
      function nIntervals()
*
      implicit none
*
      include "../commons/grid.h"
**
*     Output Variables
*
      integer nIntervals
*
      nIntervals = nin(0)
*
      return
      end
