************************************************************************
*
*     SetExternalGrid.f:
*
*     This subroutine sets the i-th x-space grid to the the user given
*     vector x.
*
************************************************************************
      subroutine SetExternalGrid(i,np,deg,x)
*
      implicit none
*
      include "../commons/grid.h"
*
*     Variables
*
      integer i,np,deg
      integer ix
      double precision x(0:nint_max)
      double precision eps
      parameter(eps=1d-12)
*
      nin(i)           = np
      inter_degree(i)  = deg
      IsExt(i)         = .true.
      ThereAreExtGrids = .true.
      do ix=0,np
         xgext(i,ix) = x(ix)
      enddo
      xmin(i) = x(0)
      if(i.eq.ngrid) xmin(i+1) = xmax + eps
*
      return
      end
