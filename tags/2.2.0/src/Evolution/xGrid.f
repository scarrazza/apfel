************************************************************************
*
*     xGrid.f:
*
*     This function returns the node values of the joint x-space grid.
*
************************************************************************
      function xGrid(alpha)
*
      implicit none
*
      include "../commons/grid.h"
**
*     Input Variables
*
      integer alpha
**
*     Output Variables
*
      double precision xGrid
*
*     Check consistency of the input variables
*
      if(alpha.lt.0.or.alpha.gt.nin(0))then
         write(6,*) "In xGrid.f:"
         write(6,*) "Invalid index, alpha =",alpha
         call exit(-10)
      endif
*
      xGrid = xg(0,alpha)
*
      return
      end
