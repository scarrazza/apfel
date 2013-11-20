************************************************************************
*
*     SetGridParameters.f:
*
*     This subroutine sets the the parameters of the i-th x-space grid.
*
************************************************************************
      subroutine SetGridParameters(i,np,deg,x)
*
      implicit none
*
      include "../commons/grid.h"
*
*     Variables
*
      integer i,np,deg
      double precision x
*
      nin(i)          = np
      inter_degree(i) = deg
      xmin(i)         = x
      if(i.eq.ngrid) xmin(i+1) = xmax
*
      return
      end
